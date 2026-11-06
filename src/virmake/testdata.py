import json
import logging
from pathlib import Path
from urllib.request import Request, urlopen

from virmake import __version__

ZENODO_RECORD_ID = "23388915"
ZENODO_API_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"

logger = logging.getLogger(__name__)


def download_test_data(destination: str | Path) -> list[Path]:
    """Download all files in the test data set from Zenodo to ``destination``."""
    logger.info("Downloading synthetic viral community...")

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)

    request = Request(
        ZENODO_API_URL,
        headers={"User-Agent": f"virmake/{__version__}"},
    )
    with urlopen(request) as response:
        record = json.load(response)

    downloaded: list[Path] = []
    for file_info in record["files"]:
        filename = Path(file_info["key"]).name
        output_path = destination / filename
        # Skip files already downloaded.
        # A size mismatch forces a re-download.
        if output_path.exists() and output_path.stat().st_size == file_info.get("size"):
            logger.info("Skipping %s (already present)", filename)
            downloaded.append(output_path)
            continue
        temporary_path = output_path.with_suffix(output_path.suffix + ".part")
        download_url = file_info["links"].get("content") or file_info["links"]["self"]
        try:
            with urlopen(download_url) as response, temporary_path.open("wb") as output:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            temporary_path.replace(output_path)
            downloaded.append(output_path)
        except Exception:
            temporary_path.unlink(missing_ok=True)
            logger.exception("Could not download %s", filename)
            raise
    return downloaded
