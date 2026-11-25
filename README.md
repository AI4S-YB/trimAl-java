# trimAl Java

**Java 11 port of trimAl v1.5.1**, with the `trimal`, `readal` and `statal` command-line
tools. The implementation is Java-only: no JNI, no invocation of the C++ binary at
runtime, and no external runtime dependencies.

The objective is to preserve the reference executable's observable behavior:
alignment bytes, statistical formatting, sequence/column selection, error messages
and exit codes. Successful differential tests establish compatibility for the tested
inputs—not an unconditional guarantee for every possible input or reference build.

## Reference version

- Upstream: https://github.com/inab/trimal
- Pinned tag: `v1.5.1`
- Commit: `d637091abe33595775f40480970d1a18d87a7bcb`
- CLI version: `trimAl v1.5.rev1 build[2025-11-25]`
- Oracle: Linux x86-64, `g++ 11.4.0`, `-Wall -O2`.

All 24 acquired C++ source/header/makefile blobs match the pinned tag. Exact source
fingerprints and the tag comparison are in `verify/reference*.json`.

This project follows the behavior-first, oracle-driven method described in the
provided tracy-java experience report. In particular, algorithms retain the
original float/double promotion points, unstable quicksort permutations and output
stream state. glibc `expf/log10f` behavior is ported explicitly rather than assumed
to match Java double arithmetic followed by a float cast.

## Build and run

Java 11+ and Maven are required for building:

```bash
mvn -f java/pom.xml package
java -jar java/target/trimal.jar --version
```

Linux, WSL or Git Bash launchers:

```bash
bash java/trimal -in alignment.fasta -out trimmed.fasta -automated1
bash java/readal -in alignment.fasta -phylip -out alignment.phy
bash java/statal -in alignment.fasta -sgc -ssc
```

Windows launchers: `java\trimal.cmd`, `java\readal.cmd`, `java\statal.cmd`.
`JAVA_HOME` is honored when set; otherwise `java` must be on `PATH`.

All three main classes are in the same JAR:

```bash
java -cp java/target/trimal.jar trimal.TrimAl -h
java -cp java/target/trimal.jar trimal.ReadAl -h
java -cp java/target/trimal.jar trimal.StatAl -h
```

Example manual trimming and selection:

```bash
bash java/trimal -in alignment.fasta -out trimmed.fasta -gt 0.9 -cons 60 -w 3
bash java/trimal -in alignment.fasta -out selected.fasta -selectcols '{' 0,2,10-15 '}'
bash java/trimal -in alignment.fasta -out trimmed.fasta -strictplus -htmlout summary.html
```

The braced selection syntax is intentional: it matches the original CLI. Options
such as `-nogap`, `-automated2`, `-colscore`, `-mw` and `-identity` do **not** exist in
this pinned reference and are rejected rather than added as compatibility aliases.

## Implemented functionality

- FASTA, CLUSTAL, NBRF/PIR, NEXUS, PHYLIP sequential/interleaved, PHYLIP-PAML and MEGA
  parsing/writing, including truncated-name variants and original whitespace.
- `-nogaps`, `-noallgaps`, `-gappyout`, `-strict`, `-strictplus`, `-automated1`.
- Manual gap/similarity/consistency thresholds, minimum coverage, half-windows,
  block sizes, terminal boundaries, complementary selection and column mappings.
- Sequence overlap filtering, representative clustering and identity thresholds.
- Multiple-alignment consistency comparisons, forced reference selection, default
  and user-defined similarity matrices.
- CDS backtranslation and original stop-codon handling options.
- Original HTML summary/color output.
- `readal` format/type/info/reverse/ungapped conversion and `statal` statistics.

## Verification

The verification scripts preserve raw stdout, stderr, exit status and **every**
produced file. No whitespace or line-ending normalization is used. Complete results
and artifact SHA-256 values are saved locally under `verify/results/` and each
independent suite's directory.

```bash
bash verify/run_all.sh
```

This command builds with unit tests enabled, runs the generated quick, extended and
legacy matrices and the companion-tool suite (`--scale` adds the slow conservation
jobs). A real test failure or timeout fails the run;
known native undefined-behavior cases are reported separately, never counted as
byte-identical PASS. See [verify/README.md](verify/README.md) for prerequisites,
reference setup, manifests and individual runners.

The integrated counts of the October 16 run are recorded in
[verify/validation-summary.json](verify/validation-summary.json); the other rows come from
their own run directories (see [Evidence provenance](#evidence-provenance)). No completed
run has an actual failure:

| Suite | Strict byte/bit matches | Separately recorded differences |
|---|---:|---:|
| Main `trimal` CLI — quick matrix | 502 / 508 | 6 native UB/crash paths |
| Main `trimal` CLI — extended matrix | 471 / 475 | 4 native UB/crash paths |
| Main `trimal` CLI — legacy baseline | 113 / 117 | 4 native UB/crash paths |
| Main `trimal` CLI — scale (3,583-seq real data + 4×60,000-col wide) | 8 / 8 | None |
| `readal` / `statal` CLI | 474 / 475 | 1 native double-free |
| Boundary/alignment audit | 480 / 605 | 125 native UB/crash paths |
| Loader fuzz (seed 90210, 180 mutations) | 180 / 180 | None (75 cases load and write; 105 exercise rejection paths) |
| Downloaded real data: Pfam-A.seed top-10 families | 40 / 40 | None |
| Downloaded real data: 8 eggNOG5 Bacteria clusters (4K-10K seqs) | 32 / 32 | None |
| Downloaded real data: Pfam-A.full mega-families (79K–366K seqs) | 30 / 30 | None |
| Companions on downloaded clusters (readal conversions / statal stats) | 64 / 64 | None |
| Thread determinism (4 workloads × 6 configs) | all byte-identical | None |
| Windows: full quick matrix (Windows Java vs Linux ref) | 498 / 502, plus 8 / 8 in a rerun of the affected cases | 4 harness-path failures in the full run (passed in the rerun) |
| Windows: full extended matrix (Windows Java vs Linux ref) | 471 / 474 | 3 native UB (same signature as Linux) |
| Windows: CLI smoke + companions | 96 / 96 | None |
| Numeric module oracle | 8,219 / 8,219 | None in the production fields |
| Independent float math oracle | 292,300 / 292,300 | None in the production fields |
| Integrated JUnit | 33 / 33 | None |

The main CLI now covers **1,100 differential cases** (quick + extended + legacy) with
zero actual failures, plus **180 loader-fuzz mutations** (seed 90210). The same quick
matrix was also run **cross-platform** with Windows Temurin Java against the Linux
reference binary: 498 of 502 cases matched in the full run, and the four remaining
cases — whose input lists the harness had not mapped to the Windows path — matched in
a rerun after the harness fix. A single full run on the current commit is pending. The
extended matrix adds full
input×output format grids,
manual-threshold/window/block/boundary/selection/cluster parameter sweeps, custom and
alternative similarity matrices, backtranslation codon variants, five-alignment
consistency sets, and edge-shape alignments (1 sequence, 1–2 columns, gapless,
50% gaps, clustering ties, duplicate names, 80-character headers, CRLF, missing
final newline, PIR-lookalike FASTA). It caught and fixed two real Java bugs
(a missing C++ `std::string::operator[]` null-character semantic and a null
sequence-mapping write on unaligned clustering input); the audit also completed
**18/18 actual writer-to-parser round trips** and three direct-method probes
(including 96 column-identity alignments at precision 9).
The suites overlap; counts are reported separately rather than added into one
inflated coverage number. Native crashes/UB are not included as successful matches.

Additional reports:

- `verify/numeric/`: real C++ module and single-precision math bit comparisons.
- `verify/alignment_audit/`: boundary cases, format round trips, direct-method probes.
- `verify/companions/`: readAl/statAl byte comparisons and Windows smoke tests.
- `java/target/surefire-reports/`: JUnit results.

### Scale validation

The four heaviest jobs — `strict`, `strictplus`, `automated1` and
`gappyout -htmlout` on a 3,583-sequence, ~7,400-column real EggNOG alignment
(O(sequence-count² × columns) conservation plus the HTML column walk) — are kept out
of the quick suite and run separately:

```bash
bash verify/run_all.sh --scale
```

All four have been completed **byte-identical** against the reference, and the Java
build is substantially faster on every one of them: strict 4.9 s vs 57.1 s,
strictplus 5.1 s vs 58.1 s, automated1 7.7 s vs 86.6 s, gappyout -htmlout 6.3 s vs
468.4 s (11–74×; totals ≈ 24 s vs ≈ 11 min; `scale_timings` in
[verify/validation-summary.json](verify/validation-summary.json), run
`scale_manifest-20261016-184943`). Three changes deliver this while keeping every
output byte identical: the O(sequence-count² × columns) statistics loops run in
parallel across rows/columns with per-unit serial accumulation order (`-Dtrimal.threads=1`
forces the plain serial loop), the HTML writers precompute each column once per
output block, and file-backed writers batch their writes instead of issuing one
flush syscall per line the way C++ `endl` does. Timings are recorded in
`scale_timings` in [verify/validation-summary.json](verify/validation-summary.json).
These jobs are excluded from the quick gate purely for runtime; a timeout still
counts as a failure.

### Mega-family scale (Pfam-A.full, up to 366,230 sequences)

The ultimate stress tier streams a 500 MB prefix window of **Pfam-A.full** (the
40.5 GB official alignment archive) and extracts six mega-families — the largest,
`7tm_1` (GPCR), holds **366,230 sequences × 3,489 columns ≈ 1.28 billion cells**,
34× larger than the biggest eggNOG cluster. Only O(n × L) methods apply here (the
O(n²) identity matrix is infeasible at this scale in *both* implementations):

| Family | Sequences × Columns | Native | Java |
|---|---|---:|---:|
| `7tm_1` — gappyout | 366,230 × 3,489 | 1,513 s | **56 s (27×)** |
| `7tm_1` — noallgaps | 366,230 × 3,489 | 1,517 s | **60 s (25×)** |
| `2_Hacid_dh_C` — gappyout | 126,849 × 1,500 | 226 s | **7 s (32×)** |
| all 30 mega cases | up to 366K seqs | **189 min total** | **5 min total (38×)** |

All 30 cases are **byte-identical** in the completed run
`pffull_manifest-20261121-114055` (`pffull_manifest.txt`, extracted by
`pffull_extract.py`). At this scale the reference spends 25 minutes and a sluggish
per-column `std::string` walk on `7tm_1` gappyout, while the port finishes in 56
seconds — and both write identical 1.3 GB outputs.

## Preserved quirks and explicit differences

### Preserved observable quirks

- Default output reproduces C++ significant-digit formatting; precision,
  fixed-point and left/right flags persist between statistics tables.
- Custom quicksort retains tie permutations relevant to cluster representatives;
  it is implemented iteratively so a near-sorted 100k-column input cannot exhaust
  the thread stack.
- PHYLIP3.2/PAML and NEXUS writer quirks, original title lines, PIR terminators,
  row grouping, trailing spaces, blank lines and LF endings are preserved.
- No automatic uppercasing of input sequence bytes or `.`/`-` conversion is added.
- Original CLI option order, mutual exclusions, misspelled messages and late error
  paths are retained, including some late errors returning status 0.
- Some numeric undefined behavior is emulated for this exact reference build:
  the all-above-threshold conservation read resolves to zero; AA use of the
  15×15 alternative matrix reads its observed contiguous linker layout. These are
  **reference-build observations**, not portable C++ guarantees.

### Safety divergences, not byte-compatible successes

1. **Non-FASTA/PIR trimming can crash the original.** Its `seqsInfo` pointer is NULL,
   but gap filtering/sequence selection/clustering dereferences it. Java uses safe
   empty auxiliary headers instead of reproducing a segmentation fault. Normal
   reference-success cases remain subject to strict comparisons.
2. **Empty-alignment `-colnumbering` prints an arbitrary native integer.** The original
   reads beyond a zero-length array; Java prints an empty `#ColumnsMap` line. This
   narrowly defined difference is separately classified with rc/stderr/files checked.
3. **Malformed CRLF compare lists can trigger a native double-free.** Java retains
   the same failed-load message but returns a defined error rather than adding a
   glibc abort diagnostic. The companion suite records this separately.
4. Zero-column slope/short-alignment/CDS edge cases avoid unsafe Java indexing where
   the source would perform undefined C++ memory accesses. The boundary audit
   records the actual native observations and Java result without counting them as
   successful byte matches.

These cases and their evidence are visible in the result reports; there is no
blanket normalization or SKIP rule hiding native-success/Java-missing-output bugs.

### Evidence provenance

The rows of the verification table were not all produced by the same code revision:

- **Re-validated on every push by CI** ([ci.yml](.github/workflows/ci.yml)): the quick,
  extended and legacy matrices, the `readal`/`statal` companions, the numeric module
  oracle, thread determinism and the seed-90210 loader fuzz.
- **Last run on 2026-10-16, before the loader fixes of commit `bd7d642` (2026-11-18)**:
  the scale gate (`scale_manifest-20261016-184943`), the boundary/alignment audit
  (`alignment_audit/results_final-review`), the independent float math oracle, and the
  Windows smoke and companion suites.
- **Last run on 2026-10-23/24, also before `bd7d642`**: the Windows quick matrix (full
  run `windows-matrix/20261023-232313` plus the 8-case rerun `20261023-232657` that
  followed the harness path fix), the Pfam-A.seed gate
  (`pfam_manifest-20261024-012205`) and the eggNOG gate
  (`eggnog_manifest-20261024-015135`).
- **Run at or after `bd7d642`**: companions on downloaded clusters
  (`readal_big_manifest-20261118-085938`, `statal_big_manifest-20261118-090822`), the
  Windows extended matrix (`windows-matrix/20261118-085948`) and the Pfam-A.full
  mega-families (`pffull_manifest-20261121-114055`).
- `verify/validation-summary.json` is the October 16 integrated snapshot; it does not
  include the fuzz, downloaded-data, thread-determinism or full Windows rows, and its
  recorded JAR predates the current build.

A single re-verification of every row on the current commit, with the commit and JAR
hashes recorded per row, is in progress; until it lands, read each row together with
its date.

## Project layout

```text
java/       Maven project, source, JUnit tests and launchers
verify/     pinned reference metadata, fixtures, differential runners and oracles
docs/       initial source-archaeology reports (with implementation corrections)
licenses/   retained third-party licenses
upstream/   local reference source/data snapshot; ignored by Git
```

## License

GPL-3.0-or-later; see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md). Adapted tracy-java
utilities retain their BSD-3-Clause notice. Ported glibc float algorithms retain
their LGPL-2.1-or-later/Sun notices. The corresponding licenses are included in the
JAR's `META-INF` directory. Distribution requires respecting the applicable source
and notice obligations; the reference executable is not a runtime dependency.
