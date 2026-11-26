"""Publish the text artifact from a FastQC archive."""

import argparse
import os
import tempfile
from pathlib import Path
from zipfile import ZipFile


def extract_fastqc_data(
    archive: str | Path,
    destination: str | Path,
) -> None:
    """Atomically publish the archive's single fastqc_data.txt member.

    Archive paths are never used as filesystem destinations. Content-aware
    validation belongs to workflow-coda's bio.fastqc verifier.
    """
    destination = Path(destination)

    with ZipFile(archive) as handle:
        members = [
            member
            for member in handle.infolist()
            if not member.is_dir()
            and member.filename.rsplit("/", 1)[-1] == "fastqc_data.txt"
        ]
        if len(members) != 1:
            raise ValueError(
                f"{archive}: expected exactly one fastqc_data.txt member, "
                f"found {len(members)}"
            )
        data = handle.read(members[0])

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as output:
            temporary = Path(output.name)
            output.write(data)

        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    """Extract the FastQC text artifact from a ZIP archive."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    extract_fastqc_data(args.archive, args.destination)


if __name__ == "__main__":
    main()
