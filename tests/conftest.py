# -*- coding: utf-8 -*-
"""Shared fixtures.

Note: this host sandbox denies python-created directories for later writes,
so tests write into the existing, writable ``data/_test/`` folder instead of
temporary directories (``data/_test`` is git-ignored).
"""
import os
import pathlib

import pytest

WORKSPACE = pathlib.Path(__file__).resolve().parent.parent
_D = WORKSPACE / "data" / "_test"


@pytest.fixture
def ws_results():
    """Existing writable folder for temporary test files."""
    _D.mkdir(parents=True, exist_ok=True)
    return str(_D)


def cleanup_results(pattern_prefix):
    for f in os.listdir(_D):
        if f.startswith(pattern_prefix):
            try:
                os.remove(os.path.join(_D, f))
            except OSError:
                pass

