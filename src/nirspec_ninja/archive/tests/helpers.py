"""Small file/response fixtures shared by offline workflow tests."""

import contextlib
import csv
import os
from pathlib import Path

import requests


NAME = "jw01473016001_02101_00001_nrs1_uncal.fits"


@contextlib.contextmanager
def working_directory(path):
    original = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def read_csv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


class Response:
    """A requests response double; supports header and midstream failures."""

    def __init__(self, content=b"abcdef", status=200, chunks=None, headers=None):
        self.status_code = status
        self.headers = headers if headers is not None else {"Content-Length": str(len(content))}
        self.url = "https://example.invalid/file"
        self.chunks = [content] if chunks is None else chunks
        self.closed = False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size):
        for chunk in self.chunks:
            if isinstance(chunk, Exception):
                raise chunk
            yield chunk

    def close(self):
        self.closed = True
