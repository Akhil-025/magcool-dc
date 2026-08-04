"""
Tests for main.py's _StreamToLogger -- the fix for stage output (from
modules that use plain print(), like core/validation.py and
core/validation_system.py) not being persisted to results/pipeline.log.
Does NOT run main() itself (full pipeline is ~6 minutes); tests the
stream class and its integration with contextlib.redirect_stdout directly.
"""
import logging
import contextlib

import main


def _make_capturing_logger():
    """A logger with an in-memory handler, isolated from main.py's own
    console/file handlers, so these tests don't touch results/pipeline.log
    or print to the real console."""
    logger = logging.getLogger("test_stream_to_logger")
    logger.handlers.clear()
    logger.propagate = False
    logger.setLevel(logging.INFO)
    records = []

    class _ListHandler(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    logger.addHandler(_ListHandler())
    return logger, records


def test_stream_to_logger_routes_complete_lines():
    logger, records = _make_capturing_logger()
    stream = main._StreamToLogger(logger)
    stream.write("hello\nworld\n")
    assert records == ["hello", "world"]


def test_stream_to_logger_buffers_partial_lines_until_newline():
    logger, records = _make_capturing_logger()
    stream = main._StreamToLogger(logger)
    stream.write("partial ")
    assert records == []
    stream.write("line\n")
    assert records == ["partial line"]


def test_stream_to_logger_flush_emits_trailing_partial_line():
    logger, records = _make_capturing_logger()
    stream = main._StreamToLogger(logger)
    stream.write("no trailing newline")
    assert records == []
    stream.flush()
    assert records == ["no trailing newline"]


def test_stream_to_logger_skips_empty_lines():
    """print() calls file.write(end) separately (default end="\\n"),
    which produces an empty-string write() between lines -- these should
    not become blank log records."""
    logger, records = _make_capturing_logger()
    stream = main._StreamToLogger(logger)
    stream.write("a\n")
    stream.write("")
    stream.write("b\n")
    assert records == ["a", "b"]


def test_redirect_stdout_captures_print_output_via_stream_to_logger():
    """End-to-end check of the actual fix: print() calls made while
    stdout is redirected to a _StreamToLogger reach the logger, matching
    how main()'s stage loop now wraps each stage."""
    logger, records = _make_capturing_logger()
    with contextlib.redirect_stdout(main._StreamToLogger(logger)):
        print("first line")
        print("second line")
    assert records == ["first line", "second line"]


def test_redirect_stdout_does_not_suppress_direct_logger_calls():
    """Stages that already call logger.info() directly (e.g. main.py's own
    run_baseline_sweep()) must be unaffected by stdout redirection -- those
    calls never touch sys.stdout in the first place."""
    logger, records = _make_capturing_logger()
    with contextlib.redirect_stdout(main._StreamToLogger(logger)):
        logger.info("direct logger call")
    assert records == ["direct logger call"]