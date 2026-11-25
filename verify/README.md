# Differential verification

The oracle is the original **trimAl v1.5.1** C++ executable, built with
`g++ -Wall -O2` on Linux x86-64. The 24 source/header/makefile blobs in the acquired
archive match tag `v1.5.1`, commit `d637091abe33595775f40480970d1a18d87a7bcb`.
See `reference.json`, `reference-tag-match.json` and `reference-source-sha256.json`.

## Requirements

- Java 11 or newer and Maven.
- Linux (or Ubuntu under WSL when launching from Windows Git Bash).
- Python 3.8+ and `g++`/`make` for rebuilding the reference oracles.
- A reference build in `verify/bin/`; these executables are ignored by Git.

## Set up the reference

From the workspace root on Linux/WSL:

```bash
git clone --branch v1.5.1 --depth 1 https://github.com/inab/trimal upstream
make -C upstream/source
mkdir -p verify/bin
cp upstream/source/{trimal,readal,statal} verify/bin/
```

If `upstream/` already contains the acquired source snapshot, do not clone over it.
The inherited baseline uses files under `upstream/dataset/`, so retain that directory.

## Run the tests

From Linux or Windows Git Bash:

```bash
bash verify/run_all.sh
```

The script builds the Java JAR with unit tests enabled, generates the deterministic
quick, extended and legacy matrices, and runs all three plus the companion suite
(`--scale` adds the slow conservation jobs). Every
native/Java pair uses the same working directory and input filename; `$OUT` and
`$HTML` expand to isolated output directories. All stdout, stderr, exit status and
produced alignment/HTML files are compared. **No whitespace, precision or line-ending
normalization is applied.** Every report records artifact byte sizes and SHA-256.

Test results go into timestamped directories under `verify/results/`; previous
results are not deleted. `report.json` is the machine-readable result. To inspect
unfinished or finished runs:

```bash
python3 verify/inspect_results.py verify/results/<run-directory>
```

To run a specific manifest, class or executable:

```bash
bash verify/compare_trimal.sh verify/quick_manifest.txt
python3 verify/differential.py verify/quick_manifest.txt --timeout 90 --jobs 2
python3 verify/differential.py verify/quick_manifest.txt --classes verify/integration/classes
```

The older `baseline_manifest.txt` is retained as the original 121-case baseline.
It contains deliberate invalid-CLI tests (for example `-nogap`, `-automated2` and
`-colscore`, none of which exist in this reference version). The generated quick
matrix additionally tests the **valid** options, including `-nogaps`, braced
`-selectcols`/`-selectseqs`, statistical combinations, windows, clustering,
backtranslation, all writers and error paths. `extended_manifest.txt` (seed 424242)
widens this further: every input×output format pair, threshold/window/block/boundary
and selection sweeps, cluster value ranges, custom/asymmetric matrices,
backtranslation codon variants, five-file consistency sets and edge-shape
alignments. `legacy_manifest.txt` reruns the original baseline under the strict
comparator, excluding only the four slow jobs that live in `scale_manifest.txt`.
All three manifests run as part of `bash verify/run_all.sh`.

## Scale tests are separate

```bash
bash verify/run_all.sh --scale
```

`scale_manifest.txt` holds the 3,583-sequence real-data conservation/strict and HTML
jobs excluded from the quick suite. The original conservation calculation is
O(sequence_count² × columns); the HTML writer also repeatedly traverses whole
columns. These jobs are not treated as passed unless the dedicated run finishes
and compares their complete outputs. A timeout is a **FAIL**, never a PASS or SKIP.

## Known native undefined behavior

The manifest syntax `SKIP:<reason> <args>` is restricted to the specified reason.
It is not an unconditional exemption:

- `orig-segv`: classified as **KNOWN_DIVERGENCE**, never PASS, only when the reference
  really crashes (SIGSEGV/SIGABRT) and neither side times out. Java may fail safely
  or produce a memory-safe result; every exit/stream/file difference is preserved.
  This is not an assertion that the outputs match.
- `orig-ub-garbage-columnsmap`: accepted only when exit status, stderr and all files
  match; the only permitted difference is the final integer on an otherwise
  matching `#ColumnsMap` line for an empty alignment. The native integer is an
  out-of-bounds memory read and changes between executions.

Native successful output versus absent Java output always remains a **FAIL**.
The reports preserve the original streams even for SKIP cases.

## Additional independent suites

- `verify/numeric/`: real C++ module oracles for float arithmetic, matrices,
  conservation, gaps, formatting, parser behavior and quicksort.
- `verify/alignment_audit/`: small boundary/alignment cases and direct-method probes.
- `verify/companions/`: `readal` and `statal` CLI comparisons.
- `java/src/test/`: JUnit tests for reusable compatibility components.

Each suite has its own runner and persisted results. Additional reruns from
Linux/WSL (all infer the workspace path from their own location):

```bash
python3 verify/numeric/run_numeric.py --label results
python3 verify/numeric/run_math.py
python3 verify/numeric/run_ub.py
python3 verify/alignment_audit/run_audit.py --round rerun
python3 verify/alignment_audit/run_probes.py rerun
python3 verify/threading_test.py
python3 verify/fuzz_suite.py
python3 verify/windows_smoke.py
python3 verify/validation_summary.py
```

**Thread determinism gate** (`threading_test.py`): for conservation/HTML workloads,
runs the Java build with `-Dtrimal.threads` = 1 (plain serial loop), 2, 7 and the
JVM default, twice for the default, and requires every rc/stderr/stdout/output file
to be byte-identical to the reference and to the other thread counts. This is the
guard behind the bit-exact parallelization claim.

**Loader fuzz gate** (`fuzz_suite.py`, seed 90210): 180 deterministic mutations of
real corpus alignments — truncations, byte flips, line surgery, CRLF injection,
control/high-byte insertion — each compared byte-for-byte through the same strict
harness. 75 of the 180 mutants are accepted by the loader and written out; the other
105 take its rejection paths, which are compared just as strictly. This gate caught
and fixed the `strcpy` NUL-truncation divergence in
`utils::readLine`. `run_numeric.py` and `run_math.py` fail on production-field
mismatches or failed processes; diagnostic comparisons to unused Java double math
are recorded separately, never mistaken for float compatibility failures. The UB
oracle preserves its native outputs and sanitizer reports without treating unsafe
C++ behavior as a guarantee. `run_audit.py` fails on actual differences, timeouts,
Java exceptions or unexecuted requested round trips; its original native crashes
remain separately classified. The Windows smoke script requires Windows `java.exe`
reachable from WSL; override `TRIMAL_WINDOWS_JAVA` or `--java-exe` when necessary.

`validation-summary.json` holds the integrated counts of the October 16 run; the
project README table lists each suite together with the date of its run (see its
"Evidence provenance" section). Generated inputs/results are local test assets, not
biological reference data.

## CI

`.github/workflows/ci.yml` automates the verification on every push: it clones
upstream tag v1.5.1, builds the reference with make, runs the JUnit-enabled Maven
build, the numeric oracle, the three differential matrices, the companion suite,
the thread determinism gate and the loader fuzz gate on an ubuntu runner. The
Windows smoke stays manual (it requires a Windows JVM host).

The workflow's final `Summarize` step currently has no effect in CI:
`validation_summary.py` depends on local audit results that CI does not produce, and
the step ignores its own exit status. In CI, each gate's result is the exit status of
its own step.

## Downloaded real-data gates (Pfam / eggNOG)

`pfam_extract.py` parses Pfam-A.seed.gz (EBI FTP, 194MB) and emits the largest
families as FASTA plus a 40-case manifest (gappyout / automated1 / strict / sgc
per family). `eggnog_stream.py` range-reads the eggNOG5 Bacteria
`2_trimmed_algs.tar` (3.4GB) over HTTP without downloading it all and extracts
mid-size clusters; the resulting 32-case manifest runs the same methods on
4,150-10,656-sequence real clusters. Both gates run the strict differential
harness and count toward the README totals:

```bash
python3 verify/pfam_extract.py            # needs verify/downloads/Pfam-A.seed.gz
python3 verify/differential.py verify/pfam_manifest.txt --timeout 3600 --jobs 2
python3 verify/eggnog_stream.py           # or: bash verify/pdl.sh <tar-url> <bytes> <out> 16 <workdir>
python3 verify/differential.py verify/eggnog_manifest.txt --timeout 3600 --jobs 1
```

Measured on the reference datasets (runs of 2026-10-24, before commit `bd7d642`):
every downloaded-data case byte-identical. On the eggNOG clusters the Java build is
faster on every case (strict up to 15x, automated1 up to 13x, sgc up to 49x). The
Pfam-A.seed families are small enough that both implementations finish in about
0.1 s per case, where JVM startup makes the Java build the slower of the two.

## CI runner environment

`.github/workflows/ci.yml` pins the runner to **ubuntu-22.04**. This matters: the
float bits of the C++ reference oracle's `expf`/`log10f` (used directly by
`statisticsConservation` for MDK and probed by the T cases) are relative to the
build environment's libm. The first CI run on `ubuntu-latest` (24.04, glibc 2.39)
produced 1,906 bit-level mismatches against the port — the same seed passes
8,219/8,219 on the documented reference environment (Ubuntu 22.04, glibc 2.35,
GCC 11, `-Wall -O2`), which is the baseline the pure-Java `expf`/`log10f` ports
target. Pinning the runner keeps the oracle environment identical to that
baseline. Actions deprecation warnings (Node.js 20) are cosmetic.
