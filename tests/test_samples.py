"""Compare against previous snakemake implementation"""

import sys
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from virmake.samples import (
    build_samples_table,
    get_contigs_samples,
    get_reads_samples,
    write_samples_table,
)

try:
    from snakemake.io import glob_wildcards

    HAVE_SNAKEMAKE = True
except ImportError:  # pragma: no cover
    HAVE_SNAKEMAKE = False

requires_snakemake = pytest.mark.skipif(
    not HAVE_SNAKEMAKE, reason="snakemake not importable as oracle"
)


reads_oracle_windows_skip = pytest.mark.skipif(
    sys.platform == "win32",
    reason=(
        "legacy reads patterns use explicit [^/]+ constraints which fail "
        "to compile in glob_wildcards on this platform; verified on linux CI"
    ),
)


# ------------------------------------------------------------------ oracle


def _legacy_get_reads_samples(path, qc=False):
    path = path.rstrip("/")

    (samples_fastq, f_fastq) = glob_wildcards(
        path + "/{sample,[^/]+}_{frac,[^/]+}.fastq"
    )
    (samples_gz, f_gz) = glob_wildcards(path + "/{sample,[^/]+}_{frac,[^/]+}.fastq.gz")

    if set(samples_fastq).intersection(set(samples_gz)):
        raise ValueError(
            "Error: compressed and uncompressed read files from the same sample."
        )

    tmp_df = pd.DataFrame(
        {
            "sample_id": samples_fastq + samples_gz,
            "p": [
                path + f"/{sample}_{frac}.fastq"
                for sample, frac in zip(samples_fastq, f_fastq, strict=False)
            ]
            + [
                path + f"/{sample}_{frac}.fastq.gz"
                for sample, frac in zip(samples_gz, f_gz, strict=False)
            ],
            "frac": f_fastq + f_gz,
        }
    )

    reads_table = pd.merge(
        tmp_df.loc[tmp_df["frac"] == "1", ["sample_id", "p"]].rename(
            columns={"p": "r1"}
        ),
        tmp_df.loc[tmp_df["frac"] == "2", ["sample_id", "p"]].rename(
            columns={"p": "r2"}
        ),
        how="inner",
        on=["sample_id"],
    )
    reads_table.sort_values(by=["sample_id"], inplace=True)

    if qc:
        reads_table = reads_table.rename(columns={"r1": "qc_r1", "r2": "qc_r2"})
    return reads_table


def _legacy_get_contigs_samples(path):
    path = path.rstrip("/")

    (samples,) = glob_wildcards(path + "/{sample}.fasta")

    samples = sorted(set(samples))
    return pd.DataFrame(
        {
            "sample_id": samples,
            "contigs": [path + f"/{sample}.fasta" for sample in samples],
        }
    )


# ----------------------------------------------------------------- helpers


def _make_tree(tmp_path: Path, files: list[str], subdir: str = "reads") -> Path:
    """Create ``<tmp_path>/<subdir>`` with the given entries.

    Entries ending in ``/`` are created as directories (to test that
    directory-like names are never matched); everything else is an empty file.
    """
    d = tmp_path / subdir
    d.mkdir()
    for name in files:
        if name.endswith("/"):
            (d / name[:-1]).mkdir()
        else:
            (d / name).write_text("")
    return d


def _normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize path-like cells for cross-platform oracle comparison."""
    out = df.copy()
    for col in out.columns:
        if col != "sample_id":
            out[col] = out[col].map(lambda v: str(Path(v)) if isinstance(v, str) else v)
    return out.reset_index(drop=True)


READS_TREES = [
    ["A_1.fastq", "A_2.fastq"],
    ["A_1.fastq.gz", "A_2.fastq.gz"],
    ["A_1.fastq", "A_2.fastq", "B_1.fastq.gz", "B_2.fastq.gz"],
    ["A_B_1.fastq", "A_B_2.fastq"],
    ["A_B_C_1.fastq", "A_B_C_2.fastq"],
    ["A.B_1.fastq", "A.B_2.fastq"],
    ["A_x1_1.fastq", "A_x1_2.fastq"],
    ["A_1.fastq", "A_2.fastq", "A_3.fastq"],
    ["A_1.fastq"],
    ["A_1.fastq", "A_2.fastq", "C_1.txt", "C.fasta", "notes.md"],
    [],
]

CONTIGS_TREES = [
    ["A.fasta", "B.fasta"],
    ["B.fasta", "A.fasta"],
    ["A_B.fasta"],
    ["A_1.fasta", "A_2.fasta"],
    ["A.fasta", "A.fastq", "A.fasta.gz", "notes.txt"],
    [],
]


# ---------------------------------------------------------- differential


class TestReadsDifferential:
    @requires_snakemake
    @reads_oracle_windows_skip
    @pytest.mark.parametrize(
        "files", READS_TREES, ids=lambda fs: ",".join(fs) or "<empty>"
    )
    def test_matches_legacy(self, tmp_path, files):
        d = _make_tree(tmp_path, files)
        ours = get_reads_samples(str(d))
        theirs = _legacy_get_reads_samples(d.as_posix())
        assert_frame_equal(_normalize(ours), _normalize(theirs), check_dtype=False)

    @requires_snakemake
    @reads_oracle_windows_skip
    @pytest.mark.parametrize("qc", [False, True], ids=["reads", "qc_reads"])
    def test_qc_rename_matches_legacy(self, tmp_path, qc):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"])
        ours = get_reads_samples(str(d), qc=qc)
        theirs = _legacy_get_reads_samples(d.as_posix(), qc=qc)
        assert_frame_equal(_normalize(ours), _normalize(theirs), check_dtype=False)

    @requires_snakemake
    @reads_oracle_windows_skip
    def test_mixed_gz_plain_raises_identically(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq.gz"])
        with pytest.raises(ValueError) as ours_exc:
            get_reads_samples(str(d))
        with pytest.raises(ValueError) as theirs_exc:
            _legacy_get_reads_samples(d.as_posix())
        assert str(ours_exc.value) == str(theirs_exc.value)


class TestContigsDifferential:
    @requires_snakemake
    @pytest.mark.parametrize(
        "files", CONTIGS_TREES, ids=lambda fs: ",".join(fs) or "<empty>"
    )
    def test_matches_legacy(self, tmp_path, files):
        d = _make_tree(tmp_path, files, subdir="contigs")
        ours = get_contigs_samples(str(d))
        theirs = _legacy_get_contigs_samples(d.as_posix())
        assert_frame_equal(_normalize(ours), _normalize(theirs), check_dtype=False)


# ---------------------------------------------------------- behavior pins


class TestReadsPins:
    def test_pairs_and_columns(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq", "B_1.fastq", "B_2.fastq"])
        t = get_reads_samples(str(d))
        assert list(t.columns) == ["sample_id", "r1", "r2"]
        assert t["sample_id"].tolist() == ["A", "B"]
        row = t.iloc[0]
        assert Path(row.r1).name == "A_1.fastq"
        assert Path(row.r2).name == "A_2.fastq"
        assert Path(row.r1).is_relative_to(Path(d))

    def test_gz_pairs(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq.gz", "A_2.fastq.gz"])
        t = get_reads_samples(str(d))
        assert t["sample_id"].tolist() == ["A"]
        assert Path(t.iloc[0].r1).name == "A_1.fastq.gz"

    def test_greedy_sample_wildcard(self, tmp_path):
        d = _make_tree(tmp_path, ["A_B_1.fastq", "A_B_2.fastq"])
        assert get_reads_samples(str(d))["sample_id"].tolist() == ["A_B"]

    def test_unpaired_reads_dropped(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq"])
        assert get_reads_samples(str(d)).empty

    def test_other_fracs_not_paired(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq", "A_3.fastq"])
        assert len(get_reads_samples(str(d))) == 1

    def test_other_fracs_participate_in_gz_check(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq.gz", "A_2.fastq.gz", "A_3.fastq"])
        with pytest.raises(ValueError, match="compressed and uncompressed"):
            get_reads_samples(str(d))

    def test_mixed_gz_plain_raises(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq.gz"])
        with pytest.raises(ValueError, match="compressed and uncompressed"):
            get_reads_samples(str(d))

    def test_hidden_files_skipped(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq", ".hidden_1.fastq"])
        assert get_reads_samples(str(d))["sample_id"].tolist() == ["A"]

    def test_directories_never_matched(self, tmp_path):
        d = _make_tree(tmp_path, ["B_1.fastq/", "B_2.fastq"])
        assert get_reads_samples(str(d)).empty

    def test_sorted_by_sample_id(self, tmp_path):
        d = _make_tree(
            tmp_path,
            [
                "C_1.fastq",
                "C_2.fastq",
                "A_1.fastq",
                "A_2.fastq",
                "B_1.fastq",
                "B_2.fastq",
            ],
        )
        assert get_reads_samples(str(d))["sample_id"].tolist() == ["A", "B", "C"]

    def test_qc_rename(self, tmp_path):
        d = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"])
        t = get_reads_samples(str(d), qc=True)
        assert list(t.columns) == ["sample_id", "qc_r1", "qc_r2"]

    def test_empty_dir_yields_empty_table(self, tmp_path):
        d = _make_tree(tmp_path, [])
        t = get_reads_samples(str(d))
        assert t.empty
        assert list(t.columns) == ["sample_id", "r1", "r2"]

    def test_missing_dir_yields_empty_table(self, tmp_path):
        t = get_reads_samples(str(tmp_path / "nope"))
        assert t.empty
        assert list(t.columns) == ["sample_id", "r1", "r2"]


class TestContigsPins:
    def test_basic_and_columns(self, tmp_path):
        d = _make_tree(tmp_path, ["A.fasta", "B.fasta"], subdir="contigs")
        t = get_contigs_samples(str(d))
        assert list(t.columns) == ["sample_id", "contigs"]
        assert t["sample_id"].tolist() == ["A", "B"]
        assert Path(t.iloc[0].contigs).name == "A.fasta"
        assert Path(t.iloc[0].contigs).is_relative_to(Path(d))

    def test_underscore_names_are_samples(self, tmp_path):
        d = _make_tree(tmp_path, ["A_B.fasta"], subdir="contigs")
        assert get_contigs_samples(str(d))["sample_id"].tolist() == ["A_B"]

    def test_end_anchored_extension(self, tmp_path):
        d = _make_tree(tmp_path, ["A.fasta.gz", "A.fastq"], subdir="contigs")
        assert get_contigs_samples(str(d)).empty

    def test_hidden_files_skipped(self, tmp_path):
        d = _make_tree(tmp_path, ["A.fasta", ".hidden.fasta"], subdir="contigs")
        assert get_contigs_samples(str(d))["sample_id"].tolist() == ["A"]

    def test_directories_never_matched(self, tmp_path):
        d = _make_tree(tmp_path, ["X.fasta/"], subdir="contigs")
        assert get_contigs_samples(str(d)).empty

    def test_sorted(self, tmp_path):
        d = _make_tree(tmp_path, ["C.fasta", "A.fasta", "B.fasta"], subdir="contigs")
        assert get_contigs_samples(str(d))["sample_id"].tolist() == ["A", "B", "C"]

    def test_empty_dir_yields_empty_table(self, tmp_path):
        d = _make_tree(tmp_path, [], subdir="contigs")
        t = get_contigs_samples(str(d))
        assert t.empty
        assert list(t.columns) == ["sample_id", "contigs"]

    def test_missing_dir_yields_empty_table(self, tmp_path):
        t = get_contigs_samples(str(tmp_path / "nope"))
        assert t.empty
        assert list(t.columns) == ["sample_id", "contigs"]


class TestBuildSamplesTable:
    def test_three_sources_merged(self, tmp_path):
        reads = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"])
        qc = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"], subdir="qc")
        contigs = _make_tree(tmp_path, ["A.fasta"], subdir="contigs")
        t = build_samples_table(
            reads=str(reads), qc_reads=str(qc), contigs=str(contigs)
        )
        assert t["sample_id"].tolist() == ["A"]
        row = t.iloc[0]
        assert Path(row.r1).name == "A_1.fastq"
        assert Path(row.qc_r2).name == "A_2.fastq"
        assert Path(row.contigs).name == "A.fasta"

    def test_disjoint_sources_outer_joined(self, tmp_path):
        reads = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"])
        qc = _make_tree(tmp_path, ["B_1.fastq", "B_2.fastq"], subdir="qc")
        t = build_samples_table(reads=str(reads), qc_reads=str(qc))
        assert set(t["sample_id"]) == {"A", "B"}
        assert t.loc[t.sample_id == "A", "qc_r1"].isna().all()
        assert t.loc[t.sample_id == "B", "r1"].isna().all()

    def test_no_sources_yields_empty_table_with_full_columns(self):
        t = build_samples_table()
        assert t.empty
        assert list(t.columns) == [
            "sample_id",
            "r1",
            "r2",
            "qc_r1",
            "qc_r2",
            "contigs",
        ]


class TestWriteSamplesTable:
    def test_writes_tsv_without_index(self, tmp_path):
        reads = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"])
        contigs = _make_tree(tmp_path, ["A.fasta"], subdir="contigs")
        out = tmp_path / "samples.tsv"
        t = write_samples_table(out, reads=str(reads), contigs=str(contigs))
        back = pd.read_csv(out, sep="\t")
        assert list(back.columns) == [
            "sample_id",
            "r1",
            "r2",
            "qc_r1",
            "qc_r2",
            "contigs",
        ]
        assert back.shape == t.shape
        assert Path(back.iloc[0].r1).name == "A_1.fastq"

    def test_gaps_roundtrip_as_nan(self, tmp_path):
        # Workflow contract: workflow_utils.get_qc_reads_loc treats both NaN
        # and "" as "no input"; to_csv writes gaps as empty cells which
        # read_csv reconstructs as NaN.
        reads = _make_tree(tmp_path, ["A_1.fastq", "A_2.fastq"])
        contigs = _make_tree(tmp_path, ["B.fasta"], subdir="contigs")
        out = tmp_path / "samples.tsv"
        write_samples_table(out, reads=str(reads), contigs=str(contigs))
        back = pd.read_csv(out, sep="\t")
        assert set(back["sample_id"]) == {"A", "B"}
        assert back.loc[back.sample_id == "A", "contigs"].isna().all()
        assert back.loc[back.sample_id == "B", "r1"].isna().all()

    def test_empty_table_writes_header_only(self, tmp_path):
        out = tmp_path / "samples.tsv"
        write_samples_table(out)
        back = pd.read_csv(out, sep="\t")
        assert back.empty
        assert list(back.columns) == [
            "sample_id",
            "r1",
            "r2",
            "qc_r1",
            "qc_r2",
            "contigs",
        ]
