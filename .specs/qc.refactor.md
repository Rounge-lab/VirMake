# QC Workflow and Input Validation Specification

## Scope

Refactor the QC stage’s sample validation, routing, provenance, retention, and
logging. Keep implementation minimal and centralized; preserve downstream
compatibility. Defer long-read processing, supplied-input content validation,
and final-report assembly.  Immutable sample sheet

Use columns: sample_id, read_type, read_layout, read_stage, r1, r2, contigs.

- read_type: short or long.
- read_layout: single or paired; mixed layouts across samples are supported.
- read_stage: raw or qc; qc means user-declared preprocessed reads.
- Read metadata is empty when no reads are supplied.
- Single-end samples require r1 and prohibit r2.
- Paired-end samples require both r1 and r2.
- contigs is optional and independent; supplied contigs may accompany reads.
- Never modify the sheet or populate it with generated paths.
- Discovery accepts _1/_2 and _R1/_R2; these suffixes declare paired reads.
- Reject missing counterparts, duplicate mates, and ambiguous discovery.
- Declare single-end inputs explicitly in the sheet.
- Explicit paths need not follow discovery naming conventions.
- Support plain and gzip-compressed FASTQ, including mixed compression within pairs.

## Validation and properties

Use one Pydantic loader shared by CLI preflight and the Snakefile. Use frozen
records with derived properties:

- sample.has_reads
- sample.is_paired
- sample.requires_trimming
- sample.has_supplied_contigs
- sample.read_type

Derive requires_trimming from has_reads and read_stage == "raw". It describes a
processing prerequisite, not whether trimming has run or is scheduled. For
paired reads, read_stage applies to both mates; mixed raw/preprocessed pairs
are unsupported. Supplied qc reads bypass trimming but still undergo FastQC.
Trimming-tool selection is configuration, not a sample property.

Validate unique IDs matching ^[A-Za-z0-9][A-Za-z0-9._-]*$, field combinations,
declared layouts, and supplied file existence. Always quote paths in shell
commands. Resolve relative paths against the sample-sheet directory. Require
reads or contigs for each sample. Recognize long reads but reject selected
read-processing stages that do not support them.

## Preflight and provenance

Perform preflight on every analysis invocation, including direct Snakemake
execution, before launching computational jobs. Validate current inputs,
selected-stage prerequisites, and reusable artifact availability. Reject
incompatible configuration or disabled prerequisites without a supplied
alternative. Future generated files need not exist; their declared producers
must be available. Allow setup commands with header-only sheets; reject empty
analysis runs. Defer supplied-input content validation and hashing. External
files can disappear after preflight; report subsequent failures clearly.

Preserve immutable input snapshots and invocation records under
rundir/provenance/. Record resolved paths, properties, sample-sheet checksum,
schema version, configuration, target, command, software versions, and
effective retention mode. Build the DAG from the same validated inputs used for
the snapshot. Check current declarations against the selected snapshot at
startup. Keep input provenance separate from workflow-coda result manifests.

## Routing and explicit dependencies

Use small, pure resolvers without copying, symlinking, or filesystem mutation.
qc_reads(sample) returns supplied preprocessed paths or generated
trimming-output paths. Schedule only prerequisites needed by selected targets
and configuration. Implement short-read single-end and paired-end QC
explicitly. Every required generated file must be a declared output and an
explicit consumer input. No downstream-required files may be produced solely as
side effects. Pass compressed inputs directly to supporting tools; declare any
required decompression.

## Verification and completion

Verify downstream-required generated artifacts through workflow-coda after each
stage. Use appropriate content verifiers; add missing verifiers rather than
substituting existence checks. QC verification includes generated QC reads,
applicable fastp JSON, and FastQC text. Supplied reads receive preflight
availability checks; their content validation is deferred. Downstream stages
depend on successful verification and the files they consume. Use strict
verification, expected-artifact declarations, and one writer per stage
manifest. QC produces results/qc.manifest.json through workflow-coda, not an
empty marker. Retain fastp JSON and FastQC fastqc_data.txt. Remove MultiQC
aggregation rules; final workflow-coda report generation is deferred. Apply
this verification contract to assembly when that stage is migrated.

## Outputs and retention

Store generated QC reads in results/reads/reads_qc/ and QC artifacts in
results/qc/. Supplied files remain at their original locations. Add --retention
keep|clean; new runs default to keep. Explicit CLI selection updates
configuration; omission uses the persisted setting. Direct Snakemake execution
uses that setting, subject to an explicit configuration override.

Classify outputs explicitly as final, reusable intermediate, or disposable
scratch. Always remove scratch after its last consumer, regardless of consumer
count. Keep mode retains reusable intermediates, including generated QC reads.
Clean mode removes designated intermediates after their last current-DAG
consumer. Preserve final results, retained QC reports, manifests, provenance,
logs, and benchmarks. Never delete supplied inputs or use --notemp to implement
keep mode. Do not make cleanable intermediates terminal targets merely to
retain them.

Keep mode supports staged execution and reuse. Clean mode is recommended for a
complete selected workflow. Later runs may regenerate removed intermediates and
refresh verification. Retained manifests describe prior verification, not
continued file availability. Clean execution does not retroactively sweep old
outputs.

## Logging and tests

Use unique job logs under rundir/logs/, including separate verification logs.
Explicitly capture tool stdout/stderr; Snakemake log: alone does not capture
output. Preserve invocation console output, command, and exit status. Stage
logs are single-writer summaries, never shared parallel append destinations.
Test mixed layouts, missing mates, suffix conventions, compression, and
configuration conflicts. Test strict verification, retention behavior, and
immutable supplied inputs and sample sheets.
