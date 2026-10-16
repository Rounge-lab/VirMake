import json
import logging
from pathlib import Path
from urllib.request import Request, urlopen


ZENODO_RECORD_ID = "23388915"
ZENODO_API_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"


def download_test_data(logger, destination):
    """Download all files in the test data set from Zenodo."""
    
    logger.info(f"\nDownloading synthetic viral community...\n")

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)

    request = Request(
        ZENODO_API_URL,
        headers={"User-Agent": "my-application/1.0"},
    )

    with urlopen(request) as response:
        record = json.load(response)

    downloaded_files = []

    for file_info in record["files"]:
        filename = Path(file_info["key"]).name
        output_path = destination / filename
        temporary_path = output_path.with_suffix(output_path.suffix + ".part")

        download_url = (
            file_info["links"].get("content")
            or file_info["links"]["self"]
        )

        try:
            with urlopen(download_url) as response:
                with temporary_path.open("wb") as output:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)

            temporary_path.replace(output_path)
            downloaded_files.append(output_path)

        except Exception:
            temporary_path.unlink(missing_ok=True)
            logger.exception("Could not download %s", filename)
            raise

    return downloaded_files


