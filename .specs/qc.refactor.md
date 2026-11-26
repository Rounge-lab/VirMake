# QC Workflow and Input Validation Specification

## Scope

Refactor sample validation, QC routing, provenance, retention, and logging.
Keep implementation minimal and centralized. Preserve downstream compatibility;
defer long-read processing, detailed input-content checks, and final-manifest
assembly.

## Immutable sample sheet

Use one row per sample with these columns: sample_id, read_type, read_stage,
r1, r2, contigs.

    - read_type is short or long; read_stage is raw or qc.
    - Both fields are empty for samples without reads.
    - r1 contains the first or only read file; r2 is optional.
    - contigs contains optional supplied contigs, independently of read inputs.
    - Allow supplied contigs alongside reads, for example for later mapping.
    - Never modify the sheet or populate it with generated output paths.
    - During directory discovery, accept _1/_2 and _R1/_R2 read suffixes.
    - Reject ambiguous discovery, duplicate mates, and incomplete declared pairs.
    - Explicit sample-sheet paths need not follow these filename conventions.
    - Support plain and gzip-compressed FASTQ inputs.

## Validation and derived properties

Implement one Pydantic loader shared by CLI preflight and the Snakefile. Use
frozen records and derive properties rather than accepting editable flags:

    - sample.has_reads
    - sample.is_paired
    - sample.is_trimmed
    - sample.has_supplied_contigs
    - sample.read_type

Validate unique safe sample IDs, field combinations, and supplied file paths.
Resolve relative paths against the sample-sheet directory. Require at least
reads or contigs for each sample. For paired inputs, read_stage applies to both
mates; mixed raw/QC pairs are unsupported. Derive is_trimmed from read_stage ==
"qc" only after validating the declared inputs. This is the user’s assertion
that all supplied reads are QC-ready, not proof of trimming. Single-end inputs
use the same declaration for their sole read file. Recognize long-read inputs,
but reject requested stages that do not support them. Keep processing-tool
selection separate from sample properties.

## Preflight on every invocation

Perform preflight whenever the pipeline runs, including direct Snakemake
execution. Validate current declarations, supplied file existence, and
selected-stage prerequisites. Check configured tool inclusion and reject
incompatible or disabled prerequisites. Check existing artifacts that the
selected workflow will reuse. Do not require future generated files to exist:
their declared producers must be available. Complete preflight before launching
computational jobs. Defer detailed FASTQ/FASTA content validation and
external-input content hashing. Preflight cannot guarantee that external files
will remain available after it completes; report later disappearance or storage
failures clearly.

## Provenance

Write an immutable input snapshot containing resolved paths and derived
properties. Record the sample-sheet checksum and snapshot schema version.
Record each invocation’s target, resolved configuration, versions, command, and
retention mode. Preserve history under rundir/provenance/. Construct the DAG
from the same validated inputs used for the snapshot. Check declarations
against the selected snapshot at workflow startup. Do not claim input-content
integrity before content checks are implemented. Keep input snapshots separate
from workflow-coda result manifests.

## Routing and explicit dependencies

Use small, pure resolvers with no copying, symlinking, or filesystem mutation.
qc_reads(sample) returns supplied QC paths or expected generated QC-read paths.
Every downstream-required generated file must be a declared producer output.
Every consumer must declare its required files as inputs. No required files may
exist solely as tool side effects. Respect selected targets and configured
stage/tool inclusion. Schedule read QC only for samples required by the
selected workflow. Implement short-read single-end and paired-end QC
explicitly. Pass compressed inputs directly to supporting tools. If
decompression is required, make it an explicit, validated
intermediate-producing rule.

## Stage verification and completion

Each stage validates all its downstream-required generated artifacts through
workflow-coda. Use a supported format verifier appropriate to each artifact’s
promised structure. Add a verifier where necessary; existence-only verification
is not content validation. The downstream stage must depend on successful
verification and the files it consumes. Retain QC-read inputs explicitly in the
QC verification contract, not only QC reports. Assembly must similarly verify
its downstream-required contigs. Use strict verification: failed required
artifacts fail the stage. Only declare artifacts expected for the selected
samples and configuration. Use one writer per stage manifest; never update one
manifest concurrently. QC completes by producing qc.manifest.json through
workflow-coda, not an empty marker. Retain fastp JSON and FastQC
fastqc_data.txt alongside the manifest. Implement QC verification now; apply
the same contract when assembly is migrated. Final report-manifest assembly
remains a separate, later step.

## Outputs and retention

Place generated QC reads under results/reads/reads_qc/. Place QC artifacts
under results/qc/; supplied inputs remain at their original locations. Add
--retention keep|clean, defaulting to keep. Always remove disposable scratch
after its last consumer. Keep mode retains reusable generated intermediates for
staged execution and reuse. Clean mode removes designated intermediates after
their current-DAG consumers finish. Preserve final results, retained QC
reports, manifests, provenance, logs, and benchmarks. Never delete
user-supplied inputs. Clean mode is recommended for a complete selected
workflow. Later invocations may regenerate removed intermediates and refresh
their verification. Retained manifests record prior verification; they do not
guarantee cleaned files still exist. Clean execution is not a retroactive sweep
of old outputs; do not use --notemp.

## Logging and tests

Store unique per-job logs under rundir/logs/, with separate verification logs.
Explicitly capture tool stdout/stderr; Snakemake log: alone does not capture
them. Capture each invocation’s console output, command, and exit status. Stage
logs are single-writer summaries, not shared parallel append destinations. Test
mixed inputs, both mate suffix conventions, gzip/plain reads, and missing
files. Test configuration conflicts, strict stage verification, and keep/clean
behavior. Test that supplied inputs and the sample sheet are never modified or
deleted.
