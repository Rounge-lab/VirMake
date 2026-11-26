"""Immutable sample-sheet declarations and shared input validation."""

import csv
import re
from pathlib import Path
from typing import Literal, Self

from pydantic import BaseModel, ConfigDict, model_validator

SCHEMA_VERSION = 1
SAMPLE_COLUMNS = (
    "sample_id",
    "read_type",
    "read_layout",
    "read_stage",
    "r1",
    "r2",
    "contigs",
)
_SAMPLE_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


class Sample(BaseModel):
    """A validated input declaration; never contains generated paths."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    sample_id: str
    read_type: Literal["short", "long"] | None = None
    read_layout: Literal["single", "paired"] | None = None
    read_stage: Literal["raw", "qc"] | None = None
    r1: Path | None = None
    r2: Path | None = None
    contigs: Path | None = None

    @property
    def has_reads(self) -> bool:
        return self.r1 is not None or self.r2 is not None

    @property
    def is_paired(self) -> bool:
        return self.has_reads and self.read_layout == "paired"

    @property
    def requires_trimming(self) -> bool:
        return self.has_reads and self.read_stage == "raw"

    @property
    def has_supplied_contigs(self) -> bool:
        return self.contigs is not None

    @model_validator(mode="after")
    def validate_declaration(self) -> Self:
        if _SAMPLE_ID_RE.fullmatch(self.sample_id) is None:
            raise ValueError("sample_id must match ^[A-Za-z0-9][A-Za-z0-9._-]*$")

        if not self.has_reads and not self.has_supplied_contigs:
            raise ValueError("each sample requires reads or contigs")

        metadata = (self.read_type, self.read_layout, self.read_stage)
        if not self.has_reads:
            if any(value is not None for value in metadata):
                raise ValueError("read metadata must be empty without reads")
            return self

        if any(value is None for value in metadata):
            raise ValueError("reads require read_type, read_layout, and read_stage")
        if self.r1 is None:
            raise ValueError("reads require r1")
        if self.read_layout == "single" and self.r2 is not None:
            raise ValueError("single-end reads prohibit r2")
        if self.read_layout == "paired" and self.r2 is None:
            raise ValueError("paired-end reads require r2")
        return self


def _resolve_input(value: str | None, directory: Path) -> Path | None:
    """Resolve an explicit input relative to the sample-sheet directory."""
    if value is None:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = directory / path
    path = path.resolve()
    if not path.is_file():
        raise ValueError(f"supplied input is not an existing file: {path}")
    return path


def load_sample_sheet(path: str | Path) -> tuple[Sample, ...]:
    """Load the seven-column TSV, returning immutable validated records.

    Header order may vary, but every schema column must occur exactly once.
    Empty cells represent absent values. Header-only sheets are valid here;
    analysis preflight must reject them.
    """
    sheet = Path(path).expanduser().resolve()
    samples = []
    seen: set[str] = set()

    with sheet.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if (
            header is None
            or len(header) != len(SAMPLE_COLUMNS)
            or set(header) != set(SAMPLE_COLUMNS)
        ):
            raise ValueError(
                f"{sheet}: expected exactly these columns: {', '.join(SAMPLE_COLUMNS)}"
            )

        for row in reader:
            line = reader.line_num
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{sheet}:{line}: wrong number of TSV fields")

            values = {
                name: value if value != "" else None for name, value in row.items()
            }
            try:
                for name in ("r1", "r2", "contigs"):
                    values[name] = _resolve_input(values[name], sheet.parent)
                sample = Sample.model_validate(values)
                if sample.sample_id in seen:
                    raise ValueError(f"duplicate sample_id: {sample.sample_id!r}")
            except ValueError as e:
                raise ValueError(f"{sheet}:{line}: {e}") from e

            seen.add(sample.sample_id)
            samples.append(sample)

    return tuple(samples)
