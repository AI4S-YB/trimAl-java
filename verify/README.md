# Differential verification

The oracle is the original **trimAl v1.5.1** C++ executable, built with
`g++ -Wall -O2` on Linux x86-64. The 24 source/header/makefile blobs in the acquired
archive match tag `v1.5.1`, commit `d637091abe33595775f40480970d1a18d87a7bcb`
(`reference.json`, `reference-tag-match.json`, `reference-source-sha256.json`;
`identify_reference.py` re-checks them). The method behind everything here is written
up in [docs/method.md](../docs/method.md).

## Requirements

- A JDK 11 or newer. Maven is not needed: the build goes through the wrapper.
- Linux, or Ubuntu under WSL when launching from Windows Git Bash.
- Python 3.8+ and `g++`/`make` for building the reference and the oracles.
- A reference build in `verify/bin/` (ignored by Git).
- For the Windows gates only: a Windows `java.exe` reachable from WSL, through
  `TRIMAL_WINDOWS_JAVA` or on PATH.
- For `numeric_junit` only: the JUnit jars of a Maven repository, through
  `MAVEN_REPOSITORY` (default `~/.m2/repository`; after one Maven build).

## Set up the reference

From the workspace root on Linux/WSL:

```bash
git clone --branch v1.5.1 --depth 1 https://github.com/inab/trimal upstream
test "$(git -C upstream rev-parse HEAD)" = d637091abe33595775f40480970d1a18d87a7bcb
make -C upstream/source
mkdir -p verify/bin
cp upstream/source/{trimal,readal,statal} verify/bin/
python3 verify/identify_reference.py
```

Clone inside WSL, or with `git -c core.autocrlf=false`: the suites read the files
under `upstream/dataset/` byte for byte, and a Windows checkout would give them
CRLF line endings. If `upstream/` already contains the acquired source snapshot, do
not clone over it; keep `upstream/dataset/`, which the manifests use.

## Run every gate

```bash
bash verify/run_all.sh
```

This builds the JAR with the unit tests enabled and then runs `verify/run_full.py`,
which sequences every gate and records one result per gate. Everything after the
build is the driver's:

```bash
python3 verify/run_full.py --list                   # gate names, needs and commands
python3 verify/run_full.py                         # every gate, in the order below
python3 verify/run_full.py --only quick,numeric    # a subset by name
python3 verify/run_full.py --skip pffull            # everything else
python3 verify/run_full.py --skip-missing           # leave out what this machine lacks
python3 verify/run_full.py --jar path/to/other.jar  # verify another build
```

A gate that needs something this machine does not have is refused with the reason:
the downloaded-data gates need their corpora under `verify/corpus/`, the `windows_*`
gates need a Windows `java.exe`, `numeric_junit` a Maven repository. With
`--skip-missing` such gates are recorded as skipped instead.

Before the driver starts, `run_all.sh` runs `verify/repo_check.py`; CI runs it
first as well. It checks the repository the way a clean checkout sees it, not the
way the developer's disk does:

- every file under `java/src/test/resources/` on disk is tracked. An unanchored
  `upstream/` in `.gitignore` once kept the golden inputs out of Git: the tests
  passed where the recording was made and failed on every clean checkout;
- the launchers, the shell scripts and the Maven wrapper carry mode 100755 in the
  index. On Windows `core.filemode` is off, so no script had ever received its
  executable bit, and `./mvnw` was refused on macOS;
- no tracked tool file names a machine of its own (a mounted drive, a home
  directory, a checked-out path);
- the workflow files are valid (`verify/workflow_check.py`). GitHub validates a
  workflow before it starts anything: an invalid file fails every run in zero
  seconds, with no job and no log, on every branch at once. A typo in an expression
  once did exactly that.

Before the first gate the driver records what it is about to verify: the commit,
whether the tree is clean, the JAR's SHA-256, the reference binaries' SHA-256 and
`--version` output, `uname`, the glibc, `g++`, `java` and Python versions, and the
system load before and after each gate. All of it ends up in
`verify/work/run-<timestamp>/gates.json`, next to each gate's log.

`verify/validation_summary.py` turns that record into `validation-summary.json`
and, with `--render`, rewrites the README tables from it. It refuses a run made in a
dirty tree, a run that merges different commits or JARs, and any run with a failure.

## The gates

| Gate | What it compares | Runner |
|---|---|---|
| `quick`, `extended`, `legacy` | `trimal` on generated and reference inputs: every option, every format, error paths | `differential.py` |
| `scale` | the slow O(n² × L) conservation and HTML jobs | `differential.py` |
| `fuzz`, `fuzz2`, `fuzz3` | `trimal` on mutated real files | `fuzz_suite.py` |
| `argfuzz`, `argfuzz2`, `argfuzz3` | `trimal` on random and accepted-shaped argument vectors | `argfuzz_suite.py` |
| `companions` | `readal` / `statal` on generated inputs | `companions/verify_companions.py` |
| `audit`, `probes` | edge shapes, writer-to-parser round trips, direct method probes | `alignment_audit/` |
| `numeric`, `math`, `ub`, `numeric_junit` | bit-level C++ module and float oracles; undefined-behaviour observations | `numeric/` |
| `threading` | 1, 2, 7 and default thread counts against each other and the reference | `threading_test.py` |
| `coverage` | every method the reference can reach and every CLI option exists in the port | `coverage_map.py` |
| `pfam`, `eggnog`, `pffull` | `trimal` on downloaded real alignments | `differential.py` |
| `readal_big`, `statal_big` | the companions on the downloaded eggNOG clusters | `differential.py` |
| `windows_*` | Windows `java.exe` against the Linux reference on the same cases | `windows_matrix.py`, `windows_smoke.py`, `companions/verify_windows_companions.py` |

Generated inputs, corpora and results are local test assets, not biological
reference data.

## The differential harness

```bash
python3 verify/differential.py verify/quick_manifest.txt --timeout 90 --jobs 2
bash verify/compare_trimal.sh verify/quick_manifest.txt     # same, also from Git Bash
python3 verify/inspect_results.py verify/results/<run-directory>
```

Every native/Java pair uses the same working directory and input names; `$OUT` and
`$HTML` expand to isolated output directories of the run, `$IN` to the repository
root. All stdout, stderr, exit status and produced alignment/HTML files are
compared. **No whitespace, precision or line-ending normalization is applied.**
Results go into timestamped directories under `verify/results/`; previous results
are not deleted. Each runner's last output line is `REPORT <path>`.

A manifest may name the CLI it exercises with a first line `# tool: readal` or
`# tool: statal`; otherwise it is `trimal`.

**A suite has to prove the reference ran.** Each report states how many cases the
reference exited 0 on (`VACUITY reference_ok=... error_path_only=...
dangling_option_args=...`). A suite on which it never does compares error messages
only; it fails unless `--allow-no-reference-success` is given. That guard found
three gates that had been recorded as fully passing while the reference was
rejecting every case.

### Known native undefined behavior

The manifest syntax `SKIP:<reason> <args>` is restricted to the specified reason.
It is not an unconditional exemption:

- `orig-segv`: classified as **KNOWN_DIVERGENCE**, never PASS, only when the
  reference really crashes (SIGSEGV, SIGABRT, SIGBUS, SIGFPE) and neither side
  times out. Java may fail safely or produce a memory-safe result; every
  exit/stream/file difference is preserved. This is not an assertion that the
  outputs match.
- `orig-ub-garbage-columnsmap`: accepted only when exit status, stderr and all files
  match; the only permitted difference is the final integer on an otherwise
  matching `#ColumnsMap` line for an empty alignment. The native integer is an
  out-of-bounds memory read and changes between executions.

Native successful output versus absent Java output always remains a **FAIL**.
The reports preserve the original streams even for SKIP cases.

## Fuzz gates

**Loader fuzz** (`fuzz_suite.py`). Three labels, 180 mutants each, every label a
fixed seed with its own corpus and manifest, so any label regenerates byte for byte:

| Label | Seed | Load-preserving mutants | Accepted by the reference |
|---|---:|---:|---:|
| `fuzz` | 90210 | 0 % | 75 of 180 |
| `fuzz2` | 90211 | 50 % | 137 of 180 |
| `fuzz3` | 90212 | 85 % | 167 of 180 |

The hard mutations are truncations, byte flips, line surgery, CRLF injection and
control/high bytes. The load-preserving ones rewrap sequences, edit headers, append
blank lines or alter spacing, and carry their kind in the file name, so a report
shows which kinds the reference accepts. Rejected mutants matter just as much: their
error messages are compared byte for byte. This gate caught and fixed the `strcpy`
NUL-truncation divergence in `utils::readLine`.

**Argument fuzz** (`argfuzz_suite.py`). Three labels of 400 argument vectors:

| Label | Seed | Vectors | Accepted by the reference | Reference crashes |
|---|---:|---|---:|---:|
| `argfuzz` | 7001 | random combinations of real options and value pools | 64 | 26 |
| `argfuzz2` | 7002 | a second random label | 71 | 24 |
| `argfuzz3` | 7003 | shapes the reference is known to accept, random values and inputs | 172 | 39 |

The reference crashes on about a tenth of the accepted-shape vectors: a trimming,
clustering or overlap method on a CLUSTAL, MEGA or PHYLIP input whose `seqsInfo` is
NULL, a selection or boundary that leaves nothing, a method that returned nothing
while `-htmlout` is set. Those vectors are listed in the generator (`NATIVE_CRASH`)
and carry the `SKIP:orig-segv` tag; the port gives a defined result or a clean error
on all of them. The first run of this gate found one real port bug (a consistency
threshold with no input threw instead of printing the reference's warning) and four
places where the port printed a stack trace where the reference crashes.

**Option sweep.** The end of the extended manifest runs every method, statistic,
format and modifier without an input alignment and on files with no and with one
sequence, and every value-taking option next to every method: 315 vectors.

## The golden CLI gate

`java/src/test/resources/golden/` holds what the reference printed and wrote for the
strictly matching cases of the quick, extended and companions suites:

```bash
python3 verify/record_golden.py --stats     # what would be recorded, and its size
python3 verify/record_golden.py             # record from the newest reports
```

`GoldenCliTest` runs the port in-process on the same arguments and requires the
same exit status, stdout, stderr and produced files. It needs neither the reference
binaries nor Linux: it is part of `mvn test`, so CI runs it on Linux, Windows and
macOS with JDK 11, 17 and 21. The working directory is `golden/root`, a copy of the
input files at their repository-relative paths (Maven sets it), so the arguments are
the very strings of the manifests and nothing is normalised. `CppOut.redirectStandard`
gives each command fresh standard streams, in the state a new process starts in.

## Other independent suites

Each can be run on its own from Linux/WSL; every runner infers the workspace path
from its own location and writes its report under an ignored directory.

```bash
python3 verify/numeric/run_numeric.py --label mine   # -> verify/numeric/runs/mine/
python3 verify/numeric/run_math.py
python3 verify/numeric/run_ub.py
python3 verify/numeric/run_junit.py                  # needs the JUnit jars
python3 verify/alignment_audit/run_audit.py --round rerun
python3 verify/alignment_audit/run_probes.py rerun
python3 verify/companions/verify_companions.py --java-jar java/target/trimal.jar
python3 verify/threading_test.py
python3 verify/coverage_map.py
python3 verify/corpus_manifest.py verify verify/corpus_manifests/eggnog.sha256
```

- **Numeric oracles.** `run_numeric.py` (8,219 probes) and `run_math.py` (292,300
  float values for `expf`/`log10f`, and 200,000 for the two places where the port
  calls double `pow()`) fail on any mismatch in a field the port uses. Comparisons to
  Java's double math that the port does not use are recorded as diagnostics only.
- **Undefined behaviour.** `run_ub.py` keeps the reference's own outputs and
  sanitizer reports for its undefined-behaviour sites without treating them as a
  guarantee.
- **Boundary audit.** `run_audit.py` fails on actual differences, timeouts, Java
  exceptions or a requested round trip that did not run; its native crashes remain
  separately classified. `run_probes.py` calls single methods of both
  implementations directly.
- **Thread determinism.** `threading_test.py` runs conservation/HTML workloads with
  `-Dtrimal.threads` = 1 (the plain serial loop), 2, 7 and the JVM default, twice for
  the default, and requires identical bytes everywhere.
- **Coverage.** `coverage_map.py` compares the port with the reference headers:
  every method the reference defines *and can reach* must exist under the same name
  (four `utils` helpers live under another name or inlined, and the script checks
  those symbols exist), and the option sets of the three CLIs must be equal. Methods
  the reference never calls, and declarations it never defines, are listed apart.
  The result is `coverage-map.json`.

## Downloaded real-data gates

The corpora are hundreds of megabytes and are not committed. What is committed is
the SHA-256 of every corpus file (`verify/corpus_manifests/`), so that a later
extraction can be compared file by file with the bytes an earlier run used. Each
extractor refuses to replace a corpus whose bytes differ from its committed list
unless `--update-manifest` is given.

```bash
# Pfam-A.seed: the 10 largest families, kept aligned
bash verify/pdl.sh https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz \
  194369799 verify/downloads/Pfam-A.seed.gz 8 verify/downloads/pfam_parts
python3 verify/pfam_extract.py                # needs verify/downloads/Pfam-A.seed.gz
python3 verify/differential.py verify/pfam_manifest.txt --timeout 3600 --jobs 2

# eggNOG5 Bacteria: eight pinned clusters and one truncated member. The stream
# reader takes the remote tar itself, or a downloaded prefix of it:
python3 verify/eggnog_stream.py                # HTTP Range requests against the tar
bash verify/pdl.sh http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_trimmed_algs.tar \
  428438264 verify/downloads/eggnog_window.bin 16 verify/downloads/eggnog_parts
python3 verify/eggnog_stream.py --from-file verify/downloads/eggnog_window.bin
python3 verify/eggnog_stream.py --verify       # corpus against the committed hashes
python3 verify/big_companions.py               # readal_big / statal_big manifests
python3 verify/differential.py verify/eggnog_manifest.txt --timeout 3600 --jobs 1
python3 verify/differential.py verify/readal_big_manifest.txt --timeout 3600 --jobs 2
python3 verify/differential.py verify/statal_big_manifest.txt --timeout 3600 --jobs 2

# Pfam-A.full: the mega-families of a prefix window
bash verify/pdl.sh https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz \
  524288000 verify/downloads/pffull_window.gz 16 verify/downloads/pffull_parts
python3 verify/pffull_extract.py               # needs verify/downloads/pffull_window.gz
python3 verify/differential.py verify/pffull_manifest.txt --timeout 3600 --jobs 1
```

`pdl.sh` is a parallel ranged downloader for servers that throttle each connection.

`pfam_extract.py` keeps the alignment: the insertion columns Pfam marks with `.` are
dropped and the `-` gaps of the model columns are kept. `eggnog_stream.py` walks the
tar headers of the archive as one stream and skips the data of members it does not
want; its eight members are pinned by name because the first prototype selected them
by accident of chunk boundaries. `COG0037.trimmed_alg.prefix.faa` is what that
prototype saved when the download window cut a member: the first 12,083,206 bytes of
COG0037, kept as a truncated-real-input case. `pffull_extract.py` streams the
families of the window one at a time and keeps the largest by sequence count.

## Windows gates

From WSL, with Windows `java.exe` reachable (`TRIMAL_WINDOWS_JAVA`, or on PATH):

```bash
python3 verify/windows_matrix.py --manifest verify/quick_manifest.txt --jobs 4
python3 verify/windows_matrix.py --manifest verify/extended_manifest.txt --jobs 4
python3 verify/windows_smoke.py
python3 verify/companions/verify_windows_companions.py
```

The Windows JVM reads the same input files through the same drive and writes into
its own per-case directory; everything is then compared with what the Linux
reference produced on the same arguments.

## CI

`.github/workflows/ci.yml` has two jobs.

`verify` runs on **ubuntu-22.04** on every push and pull request: it clones the
reference at the pinned commit, re-checks its 24 blobs, builds and unit-tests the
port through the Maven wrapper, regenerates the manifests and fails if they differ
from the committed ones, then runs the coverage check, both float oracles, the three
differential matrices, the companions, the six fuzz labels, the boundary audit with
its probes and the thread determinism gate. It fails if any runner modified a
tracked file, and uploads every report as an artifact.

`unit` runs the same `mvn package` on ubuntu-22.04, windows-latest and macos-latest
with JDK 11, 17 and 21. Its golden CLI gate needs no reference binary, which is what
makes the port testable on platforms the reference cannot be built for.

`.github/workflows/scale.yml` runs the scale gate weekly and on demand. The
downloaded-data and Windows gates stay local: they need corpora and a host CI does
not have.

## CI runner environment

The `verify` job is pinned to **ubuntu-22.04**. This matters: the float bits of the
C++ reference oracle's `expf`/`log10f` (used directly by `statisticsConservation` for
MDK and probed by the T cases) are relative to the build environment's libm. The
first CI run on `ubuntu-latest` (24.04, glibc 2.39) produced 1,906 bit-level
mismatches against the port — the same seed passes 8,219/8,219 on the documented
reference environment (Ubuntu 22.04, glibc 2.35, GCC 11, `-Wall -O2`), which is the
baseline the pure-Java `expf`/`log10f` ports target. Pinning the runner keeps the
oracle environment identical to that baseline.
