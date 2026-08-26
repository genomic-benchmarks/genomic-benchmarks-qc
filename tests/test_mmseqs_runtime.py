"""Tests for talking to the MMseqs2 binary.

Two halves. The preflight exists to turn one specific crash into a sentence: a
binary built for a newer CPU dies with SIGILL, which tells the user nothing.
What it must not do is stop a machine MMseqs2 would have run on perfectly well,
so the branches that give up on checking are as much the subject here as the
one that refuses.

The rest is the subprocess itself, and it is here for the failures a user
cannot diagnose from the outside - a search that exits non-zero, a search that
has to be killed, a progress bar that would otherwise fill the log file. None
of it needs MMseqs2: a fake binary on PATH can exit any way the real one can,
and can be made to hang, which the real one cannot on demand.
"""

import json
import logging
import os
import subprocess
import sys
import textwrap

import pytest

from genomic_benchmarks_qc.utils import mmseqs_runtime


@pytest.fixture
def mmseqs_on_path(monkeypatch):
    """Pretend an `mmseqs` binary is installed, so the CPU branches are reached."""
    monkeypatch.setattr(mmseqs_runtime.shutil, 'which', lambda _: '/usr/bin/mmseqs')


class TestPreflight:

    def test_a_missing_binary_is_refused(self, monkeypatch):
        monkeypatch.setattr(mmseqs_runtime.shutil, 'which', lambda _: None)
        with pytest.raises(RuntimeError, match='not found in PATH'):
            mmseqs_runtime.check_mmseqs_preflight()

    def test_a_cpu_without_any_supported_flag_is_refused(self, mmseqs_on_path, monkeypatch):
        monkeypatch.setattr(mmseqs_runtime.platform, 'system', lambda: 'Linux')
        monkeypatch.setattr(mmseqs_runtime.platform, 'machine', lambda: 'x86_64')
        monkeypatch.setattr(mmseqs_runtime, '_read_linux_cpu_flags', lambda: {'fpu', 'mmx'})
        with pytest.raises(RuntimeError, match='instruction set flags'):
            mmseqs_runtime.check_mmseqs_preflight()

    def test_an_x86_cpu_with_a_supported_flag_passes(self, mmseqs_on_path, monkeypatch):
        monkeypatch.setattr(mmseqs_runtime.platform, 'system', lambda: 'Linux')
        monkeypatch.setattr(mmseqs_runtime.platform, 'machine', lambda: 'x86_64')
        monkeypatch.setattr(mmseqs_runtime, '_read_linux_cpu_flags', lambda: {'avx2', 'sse2'})
        mmseqs_runtime.check_mmseqs_preflight()

    def test_a_non_linux_system_is_let_through(self, mmseqs_on_path, monkeypatch, caplog):
        monkeypatch.setattr(mmseqs_runtime.platform, 'system', lambda: 'Darwin')
        with caplog.at_level(logging.WARNING):
            mmseqs_runtime.check_mmseqs_preflight()
        assert 'Darwin' in caplog.text

    @pytest.mark.parametrize('arch', ['aarch64', 'arm64', 'ppc64le'])
    def test_a_non_x86_linux_is_let_through_with_a_warning(
            self, mmseqs_on_path, monkeypatch, caplog, arch):
        """MMseqs2 ships ARM builds, so refusing here refuses working hardware.

        The flags the check looks for are x86_64 ones. On an architecture that
        has none of them there is nothing to check, and nothing to check is not
        the same as something wrong.
        """
        monkeypatch.setattr(mmseqs_runtime.platform, 'system', lambda: 'Linux')
        monkeypatch.setattr(mmseqs_runtime.platform, 'machine', lambda: arch)

        def unreadable():
            raise AssertionError('the x86 CPU flags must not be consulted on ' + arch)

        monkeypatch.setattr(mmseqs_runtime, '_read_linux_cpu_flags', unreadable)

        with caplog.at_level(logging.WARNING):
            mmseqs_runtime.check_mmseqs_preflight()

        assert arch in caplog.text

    def test_unreadable_cpu_flags_on_x86_are_still_refused(
            self, mmseqs_on_path, monkeypatch):
        """The one architecture the check does cover has to stay covered."""
        monkeypatch.setattr(mmseqs_runtime.platform, 'system', lambda: 'Linux')
        monkeypatch.setattr(mmseqs_runtime.platform, 'machine', lambda: 'AMD64')
        monkeypatch.setattr(mmseqs_runtime, '_read_linux_cpu_flags', lambda: None)
        with pytest.raises(RuntimeError, match='/proc/cpuinfo'):
            mmseqs_runtime.check_mmseqs_preflight()


# What every fake gets before its body: the command line it was called with, and
# the modules a body might want.
FAKE_PREAMBLE = '''
import json, os, sys, time
argv = sys.argv[1:]
'''


def fake_mmseqs(tmp_path, monkeypatch, body):
    """Put an `mmseqs` on PATH that runs `body`, and return the directory it is in.

    `body` is Python, run with this interpreter, with `argv` bound to the
    command line MMseqs2 was called with - so `argv[3]` is the output table it
    is expected to write. Nothing about the tests needs the real binary; what
    they are about is how `run_search` behaves when it exits, or does not.
    """
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    script = bin_dir / 'mmseqs'
    script.write_text(f'#!{sys.executable}\n' + FAKE_PREAMBLE + textwrap.dedent(body))
    script.chmod(0o755)
    monkeypatch.setenv('PATH', str(bin_dir) + os.pathsep + os.environ['PATH'])
    return bin_dir


# Enough of a search to succeed: record the call, say something, write the table.
SUCCEEDS = r"""
    open(os.environ['CALL_RECORD'], 'w').write(
        json.dumps({'argv': argv, 'tty': os.environ.get('TTY')}))
    print('easy-search: 3 hits')
    open(argv[3], 'w').write('query\ttarget\n')
"""


@pytest.fixture
def call_record(tmp_path, monkeypatch):
    """Where the fake writes what it was called with."""
    path = tmp_path / 'call.json'
    monkeypatch.setenv('CALL_RECORD', str(path))
    return path


@pytest.fixture
def search_args(tmp_path):
    """The four paths every `run_search` call needs, ready to splat."""
    (tmp_path / 'test.fa').write_text('>seq_0_test\nACGT\n')
    (tmp_path / 'train.fa').write_text('>seq_0_train\nACGT\n')
    return {
        'test_fasta_file': tmp_path / 'test.fa',
        'train_fasta_file': tmp_path / 'train.fa',
        'out_file': tmp_path / 'hits.tsv',
        'tmp_dir': tmp_path / 'scratch',
    }


class TestASearchThatWorks:

    def test_it_returns_the_table_mmseqs_wrote(self, tmp_path, monkeypatch, search_args,
                                               call_record):
        fake_mmseqs(tmp_path, monkeypatch, SUCCEEDS)

        result = mmseqs_runtime.run_search(**search_args)

        assert result == search_args['out_file']
        assert result.read_text() == 'query\ttarget\n'

    def test_the_command_line_is_one_the_summariser_can_read_back(
            self, tmp_path, monkeypatch, search_args, call_record):
        """`--format-mode 4` is the header row the chunked reader joins on."""
        fake_mmseqs(tmp_path, monkeypatch, SUCCEEDS)

        mmseqs_runtime.run_search(**search_args)
        argv = json.loads(call_record.read_text())['argv']

        assert argv[0] == 'easy-search'
        assert argv[1:5] == [str(search_args['test_fasta_file']),
                             str(search_args['train_fasta_file']),
                             str(search_args['out_file']),
                             str(search_args['tmp_dir'])]
        assert argv[argv.index('--format-output') + 1] == ','.join(
            mmseqs_runtime.MMSEQS_REQUIRED_COLS)
        assert argv[argv.index('--format-mode') + 1] == '4'
        assert argv[argv.index('--search-type') + 1] == '3'
        assert argv[argv.index('--strand') + 1] == '1'

    def test_the_tuning_options_are_passed_only_when_they_are_given(
            self, tmp_path, monkeypatch, search_args, call_record):
        fake_mmseqs(tmp_path, monkeypatch, SUCCEEDS)

        mmseqs_runtime.run_search(**search_args)
        assert '--threads' not in json.loads(call_record.read_text())['argv']

        mmseqs_runtime.run_search(**search_args, threads=4, split_memory_limit='8G')
        argv = json.loads(call_record.read_text())['argv']
        assert argv[argv.index('--threads') + 1] == '4'
        assert argv[argv.index('--split-memory-limit') + 1] == '8G'

    def test_mmseqs_is_told_it_has_a_terminal(self, tmp_path, monkeypatch, search_args,
                                              call_record):
        """Without TTY=1 it prints no progress, and a long search goes silent."""
        fake_mmseqs(tmp_path, monkeypatch, SUCCEEDS)

        mmseqs_runtime.run_search(**search_args)

        assert json.loads(call_record.read_text())['tty'] == '1'


class TestWhatReachesTheLog:

    def test_a_progress_bar_is_logged_once_per_interval(self, tmp_path, monkeypatch,
                                                        search_args, caplog):
        """MMseqs2 redraws one line thousands of times over a long search."""
        fake_mmseqs(tmp_path, monkeypatch, r"""
            for i in range(200):
                sys.stdout.write('[=====] %d%%\r\n' % i)
            sys.stdout.flush()
            open(argv[3], 'w').write('')
        """)

        with caplog.at_level(logging.DEBUG, logger='genomic_benchmarks_qc'):
            mmseqs_runtime.run_search(**search_args)

        kept = [r for r in caplog.records if 'MMSeqs2 progress' in r.getMessage()]
        assert len(kept) == 1, f'{len(kept)} of 200 redraws reached the log'

    def test_ordinary_output_is_not_throttled(self, tmp_path, monkeypatch, search_args,
                                              caplog):
        fake_mmseqs(tmp_path, monkeypatch, r"""
            for i in range(5):
                print('Stage %d complete' % i)
            open(argv[3], 'w').write('')
        """)

        with caplog.at_level(logging.DEBUG, logger='genomic_benchmarks_qc'):
            mmseqs_runtime.run_search(**search_args)

        assert sum('Stage' in r.getMessage() for r in caplog.records) == 5

    def test_what_mmseqs_says_on_stderr_is_logged_as_an_error(self, tmp_path, monkeypatch,
                                                             search_args, caplog):
        """The reason is on stderr; the exception carries none of it."""
        fake_mmseqs(tmp_path, monkeypatch, r"""
            sys.stderr.write('Invalid database name\n')
            sys.exit(1)
        """)

        with caplog.at_level(logging.ERROR, logger='genomic_benchmarks_qc'), \
                pytest.raises(RuntimeError, match='MMSeqs2 search failed'):
            mmseqs_runtime.run_search(**search_args)

        assert 'Invalid database name' in caplog.text
        assert 'return code: 1' in caplog.text


class TestASearchThatDoesNotFinish:

    def test_a_non_zero_exit_is_raised_rather_than_returned(self, tmp_path, monkeypatch,
                                                            search_args):
        """A failed search can still leave an output file, so the code is the answer."""
        fake_mmseqs(tmp_path, monkeypatch, r"""
            open(argv[3], 'w').write('')
            sys.exit(3)
        """)

        with pytest.raises(RuntimeError, match='MMSeqs2 search failed'):
            mmseqs_runtime.run_search(**search_args)

    def test_an_interrupted_run_does_not_leave_mmseqs_running(self, tmp_path, monkeypatch,
                                                              search_args):
        """A search is the whole machine for as long as it runs.

        Ctrl-C arrives as KeyboardInterrupt, which is not an Exception, so
        catching Exception here left the search running in the case where it
        cost the most. The fake hangs, so there is really a process to kill.
        """
        fake_mmseqs(tmp_path, monkeypatch, r"""
            time.sleep(120)
        """)

        started, interrupted = [], []

        class InterruptingPopen(subprocess.Popen):
            """A Popen whose first `wait` is Ctrl-C, and whose second is real."""

            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                started.append(self)

            def wait(self, *args, **kwargs):
                if not interrupted:
                    interrupted.append(True)
                    raise KeyboardInterrupt('ctrl-c')
                return super().wait(*args, **kwargs)

        monkeypatch.setattr(mmseqs_runtime.subprocess, 'Popen', InterruptingPopen)

        with pytest.raises(KeyboardInterrupt):
            mmseqs_runtime.run_search(**search_args)

        assert len(started) == 1
        assert started[0].poll() is not None, 'mmseqs is still running'

    def test_a_binary_that_is_not_there_stops_the_run_before_it_starts(
            self, tmp_path, monkeypatch, search_args):
        (tmp_path / 'empty').mkdir()
        monkeypatch.setenv('PATH', str(tmp_path / 'empty'))

        with pytest.raises(RuntimeError, match='not found in PATH'):
            mmseqs_runtime.run_search(**search_args)

        assert not search_args['out_file'].exists()
