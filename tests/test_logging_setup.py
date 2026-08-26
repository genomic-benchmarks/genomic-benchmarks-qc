"""Tests for whose logging this package is allowed to configure.

A library that calls `logging.basicConfig` configures the root logger, which
belongs to whoever imported it: every other library in that process starts
reporting through handlers it never asked for. `basicConfig` also does nothing
once the root logger has handlers, which made the second `run()` in a process
ignore its own `log_level` and `log_file` and keep writing to the first one's.

Everything here logs to `genomic_benchmarks_qc`, and `setup_logger` configures
that and nothing else.
"""

import logging
import subprocess
import sys

import pytest

from genomic_benchmarks_qc.utils.input_utils import PACKAGE_LOGGER_NAME, setup_logger


@pytest.fixture
def package_logger():
    """The package logger, restored afterwards - `setup_logger` is global state."""
    package_logger = logging.getLogger(PACKAGE_LOGGER_NAME)
    handlers, level = package_logger.handlers[:], package_logger.level
    yield package_logger
    for handler in package_logger.handlers[:]:
        package_logger.removeHandler(handler)
    for handler in handlers:
        package_logger.addHandler(handler)
    package_logger.setLevel(level)


@pytest.fixture
def root_logger():
    """The caller's root logger, watched for anything the package does to it."""
    root = logging.getLogger()
    handlers, level = root.handlers[:], root.level
    yield root
    root.handlers[:] = handlers
    root.setLevel(level)


class TestTheRootLoggerIsLeftAlone:

    def test_importing_the_package_configures_nothing(self):
        """In a subprocess: by now this process has imported it many times over."""
        script = (
            "import logging\n"
            "root = logging.getLogger()\n"
            "before = (len(root.handlers), root.level)\n"
            "import genomic_benchmarks_qc.evaluate_classes\n"
            "print((len(root.handlers), root.level) == before)\n"
        )
        result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)

        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == 'True'

    def test_setting_the_log_up_adds_no_handler_to_it(self, package_logger, root_logger,
                                                      tmp_path):
        before = list(root_logger.handlers)

        setup_logger(logging.DEBUG, tmp_path / 'run.log')

        assert list(root_logger.handlers) == before

    def test_the_matplotlib_logger_is_not_pinned_either(self, package_logger, root_logger):
        """It was, and only because a root handler at DEBUG made it necessary."""
        matplotlib_logger = logging.getLogger('matplotlib')
        before = matplotlib_logger.level

        setup_logger(logging.DEBUG)

        assert matplotlib_logger.level == before

    def test_an_application_watching_the_root_still_sees_the_records(
            self, package_logger, root_logger, caplog):
        """Leaving the root alone is not the same as talking past it."""
        setup_logger(logging.INFO)

        with caplog.at_level(logging.INFO):
            logging.getLogger(PACKAGE_LOGGER_NAME + '.somewhere').info('a record')

        assert 'a record' in caplog.text


class TestSettingItUpTwice:

    def test_the_second_level_is_the_one_that_applies(self, package_logger, tmp_path):
        setup_logger(logging.DEBUG, tmp_path / 'first.log')
        setup_logger(logging.WARNING, tmp_path / 'second.log')

        logging.getLogger(PACKAGE_LOGGER_NAME).info('an INFO record')

        assert 'an INFO record' not in (tmp_path / 'second.log').read_text()

    def test_the_second_file_is_the_one_written_to(self, package_logger, tmp_path):
        setup_logger(logging.INFO, tmp_path / 'first.log')
        logging.getLogger(PACKAGE_LOGGER_NAME).info('from the first run')
        setup_logger(logging.INFO, tmp_path / 'second.log')
        logging.getLogger(PACKAGE_LOGGER_NAME).info('from the second run')

        first = (tmp_path / 'first.log').read_text()
        second = (tmp_path / 'second.log').read_text()
        assert 'from the first run' in first and 'from the second run' not in first
        assert 'from the second run' in second and 'from the first run' not in second

    def test_the_handlers_are_replaced_rather_than_piled_up(self, package_logger, tmp_path):
        """Otherwise each run's lines appear once per run that came before it."""
        for index in range(4):
            setup_logger(logging.INFO, tmp_path / f'{index}.log')

        assert len(package_logger.handlers) == 2

    def test_a_run_without_a_file_stops_writing_to_the_last_one(self, package_logger, tmp_path):
        setup_logger(logging.INFO, tmp_path / 'run.log')
        setup_logger(logging.INFO)

        logging.getLogger(PACKAGE_LOGGER_NAME).info('after the file was dropped')

        assert 'after the file was dropped' not in (tmp_path / 'run.log').read_text()
