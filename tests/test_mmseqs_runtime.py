"""Tests for the MMseqs2 preflight check.

The preflight exists to turn one specific crash into a sentence: a binary built
for a newer CPU dies with SIGILL, which tells the user nothing. What it must not
do is stop a machine MMseqs2 would have run on perfectly well, so the branches
that give up on checking are as much the subject here as the one that refuses.
"""

import logging

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
