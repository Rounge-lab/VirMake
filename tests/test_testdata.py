"""Tests for the test-data downloader."""

import json

from virmake.testdata import download_test_data


class _FakeResponse:
    def __init__(self, payload: bytes):
        self._payload = payload
        self._consumed = False

    def read(self, size: int = -1) -> bytes:
        if self._consumed:
            return b""
        self._consumed = True
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _record(files: list[dict]) -> bytes:
    return json.dumps({"files": files}).encode()


def test_skips_existing_file_with_matching_size(tmp_path, monkeypatch):
    existing = tmp_path / "synthetic_1.fastq.gz"
    existing.write_bytes(b"x" * 100)
    calls = []

    def fake_urlopen(request, **kwargs):
        calls.append(request)
        return _FakeResponse(
            _record(
                [
                    {
                        "key": "synthetic_1.fastq.gz",
                        "size": 100,
                        "links": {"content": "https://example.com/f"},
                    }
                ]
            )
        )

    monkeypatch.setattr("virmake.testdata.urlopen", fake_urlopen)
    assert download_test_data(tmp_path) == [existing]
    assert len(calls) == 1  # API call only — no download attempt
    assert existing.read_bytes() == b"x" * 100  # content untouched


def test_size_mismatch_forces_redownload(tmp_path, monkeypatch):
    existing = tmp_path / "synthetic_1.fastq.gz"
    existing.write_bytes(b"x")  # size 1, record says 2
    calls = []

    def fake_urlopen(request, **kwargs):
        calls.append(request)
        if len(calls) == 1:
            return _FakeResponse(
                _record(
                    [
                        {
                            "key": "synthetic_1.fastq.gz",
                            "size": 2,
                            "links": {"content": "https://example.com/f"},
                        }
                    ]
                )
            )
        return _FakeResponse(b"zz")

    monkeypatch.setattr("virmake.testdata.urlopen", fake_urlopen)
    download_test_data(tmp_path)
    assert existing.read_bytes() == b"zz"
