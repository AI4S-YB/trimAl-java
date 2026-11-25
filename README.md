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
- Oracle: Linux x86-64, `g++ 11.4.0`, `-Wall -O2` on Ubuntu 22.04.

All 24 acquired C++ source/header/makefile blobs match the pinned tag. Exact source
fingerprints and the tag comparison are in `verify/reference*.json`, and CI re-checks
the blobs on every run.

The reference's float bits are relative to the environment it was built in, and that
was measured ([verify/reference-environments.json](verify/reference-environments.json)):
a reference built on Ubuntu 24.04 (glibc 2.39, GCC 13) prints the same bits as one
built on 22.04, and the port matches both, probe for probe. On Ubuntu 26.04 (glibc
2.43, GCC 15) it does not: `log10f` became a correctly rounded implementation, and
GCC 15 contracts one accumulation into a fused multiply-add. The port is the 22.04
reference, which is also the 24.04 one, and CI is pinned to that image.

The port is built and verified by a behavior-first, oracle-driven method, written up
in [docs/method.md](docs/method.md): algorithms retain the original float/double
promotion points, unstable quicksort permutations and output stream state, and glibc
`expf/log10f` behavior is ported explicitly rather than assumed to match Java double
arithmetic followed by a float cast.

## Build and run

A JDK 11 or newer is required; Maven is not, the build goes through the Maven
wrapper:

```bash
java/mvnw -f java/pom.xml package        # Windows: java\mvnw.cmd -f java/pom.xml package
java -jar java/target/trimal.jar --version
```

`package` also runs the unit tests, a golden gate of about 1,700 recorded reference
cases among them; skip them with `-DskipTests`. Two builds of one commit produce the
same JAR bytes.

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

Equality here means equality of everything a run produces: exit status, stdout,
stderr and every written file, compared byte for byte. No whitespace, precision or
line-ending normalization is applied anywhere.

### One driver, one record

```bash
bash verify/run_all.sh
```

builds the JAR with the unit tests enabled and runs `verify/run_full.py`, which
sequences every gate on that one build. Before the first gate the driver records the
commit, whether the tree was clean, the JAR's SHA-256, the reference binaries'
SHA-256 and versions and the toolchain; per gate it keeps the log, the exit status,
the seconds and the system load. `verify/validation_summary.py` turns that record into
[verify/validation-summary.json](verify/validation-summary.json) and renders the
tables below from it, so a number in them is never typed by hand. The summary refuses
a run made in a dirty tree, runs that verified different JARs or different reference
binaries, and any run with a failure. It does accept runs of several commits when the
build is the same — after a change to a corpus, say, only the gates that read it are
re-run — and then every gate names the commit and the run it comes from.
[verify/README.md](verify/README.md) lists every gate and shows how to run one on its
own.

### Results

<!-- verification-table:start -->
| Suite | Strict byte/bit matches | Separately recorded differences | Reference accepted / rejected |
|---|---:|---|---|
| Main `trimal` CLI — quick matrix | 502 / 508 | 6 native UB/crash paths | 441 / 67 |
| Main `trimal` CLI — extended matrix | 779 / 790 | 11 native UB/crash paths | 657 / 133 |
| Main `trimal` CLI — legacy baseline | 113 / 117 | 4 native UB/crash paths | 76 / 41 |
| Main `trimal` CLI — scale (3,583-sequence real data + 4×60,000-column wide) | 8 / 8 † | None | 8 / 0 |
| `readal` / `statal` CLI | 474 / 475 | 1 native UB/crash path | — |
| Boundary/alignment audit (parse, writers, edge shapes, round trips) | 480 / 605 | 125 native UB/crash paths | — |
| Direct-source probes (matrix order, locale, 96 identity alignments) | 3 / 3 | None | — |
| Loader fuzz (three labels × 180 mutations) | 540 / 540 | None | 379 / 161 |
| CLI argument fuzz (two random labels and one accepted-shape label, 400 vectors each) | 1,111 / 1,200 | 89 native UB/crash paths | 307 / 804 |
| Downloaded real data: Pfam-A.seed top-10 families | pending † | — | — |
| Downloaded real data: 8 eggNOG5 Bacteria clusters + 1 truncated member | pending † | — | — |
| Downloaded real data: Pfam-A.full mega-families (up to 366K sequences) | 30 / 30 † | None | 30 / 0 |
| Companions on downloaded clusters (readal conversions / statal stats) | 90 / 90 † (readal only) | None | 84 / 6 |
| Thread determinism (4 workloads × 6 thread settings) | all byte-identical | None | — |
| Windows: full quick matrix (Windows Java vs Linux reference) | 499 / 502 in the full run, plus 3 / 3 rerun in isolation † | 3 load-induced failures in the full run (passed in isolation) | — |
| Windows: full extended matrix (Windows Java vs Linux reference) | 471 / 474 † | 3 native UB/crash paths | — |
| Windows: CLI smoke + companions (Windows Java vs Linux reference) | 96 / 96 † | None | — |
| Numeric module oracle (bit-level, incl. parallel paths) | 8,219 / 8,219 | None | — |
| Float math oracle (bit-level: the float overloads and the two double pow() sites) | 492,300 / 492,300 | None | — |
| Native undefined-behaviour observations (repeated native runs vs Java) | 18 / 18 | 15 native UB/crash paths | — |
| Reference coverage check (reachable methods and CLI options) | 140 / 140 | None | — |
| Golden CLI gate: recorded reference output, port run in-process in `mvn test` | 1,432 / 1,432 | None | — |
| Integrated JUnit (`mvn package`), the golden gate included | 1,466 / 1,466 | None | — |
<!-- verification-table:end -->

<!-- verification-provenance:start -->
**Provenance (provisional).** Rows marked † are the last completed run of that gate
before this build: the scale gate (2026-10-16), the Pfam-A.full gate (2026-11-21),
the Windows suites (2026-10-23 and 2026-11-18) and the `readal` half of the
companions on the eggNOG clusters (2026-11-25, before the last fixes). The Pfam-A.seed
and eggNOG gates had not yet been run on a corpus the reference accepts. A full run
of every gate on one commit is in progress; this paragraph is rendered from its
record when it has finished.
<!-- verification-provenance:end -->

### How to read the table

- A *strict match* is a case whose reference run and port run produced the same exit
  status, the same stdout and stderr and the same files, byte for byte.
- *Native UB/crash paths* are cases in which the reference reads past an array or
  crashes. What the port does there is a defined behaviour, listed under
  [safety divergences](#safety-divergences-not-byte-compatible-successes); it is not
  a match and is never counted as one. A vector carries its crash tag only for as
  long as the reference really crashes on it: on a run where it exits normally, the
  vector is compared strictly.
- For the fuzz suites the table also states how many vectors the reference accepted.
  Every differential report carries that count. A suite on which the reference never
  exits 0 would compare error messages only, so it fails unless that is declared on
  the command line. That guard found three gates that had been recorded as passing
  while the reference rejected every case: the Pfam-A.seed corpus had been written
  without its gaps, so trimAl saw unaligned sequences, and two companion manifests
  used options the tools do not have.
- The suites overlap; the counts are reported per suite and never summed into one
  coverage figure.
- ‡ Six cases of the boundary audit — a backtranslation whose selection of columns or
  sequences is empty — sit on undefined behaviour on which the reference crashes in
  about half of the runs and exits normally in the others (measured: 10 to 13 crashes
  in 20 runs, per case). Whenever it exits normally, its output is byte-identical to
  the port's. Those cases are therefore strict matches in some runs and native crash
  paths in others; the two together are always 605.

### What the gates are

- **Differential matrices.** `trimal` on generated and reference inputs. The quick
  matrix keeps the original error-path baseline and adds every valid option; the
  extended matrix adds full input×output format grids, parameter sweeps, custom and
  alternative similarity matrices, backtranslation codon variants, five-alignment
  consistency sets and edge-shape alignments (1 sequence, 1–2 columns, gapless, 50 %
  gaps, clustering ties, duplicate names, 80-character headers, CRLF, missing final
  newline), and finally every method, statistic, format and modifier without an input
  alignment and on files with no and with one sequence.
- **Fuzz.** Three labels of 180 mutated real files (truncations, byte flips, line
  surgery, CRLF, control bytes, and load-preserving edits such as rewrapped sequences
  and edited headers) and three labels of 400 argument vectors. Rejected inputs
  matter as much as accepted ones: their error messages are compared byte for byte.
- **Golden CLI gate.** The reference's own output for the strictly matching cases is
  recorded under `java/src/test/resources/golden/`, and `mvn test` runs the port
  in-process on the same arguments. It needs neither the reference binaries nor
  Linux, which is what makes the port testable on platforms the reference cannot be
  built for; CI runs it on Linux, Windows and macOS with JDK 11, 17 and 21.
- **Oracles.** Bit-level comparisons against the real C++ modules, the float
  overloads of `expf`/`log10f` over 292,300 values, and the two places where the
  port calls double `pow()` over 200,000.
- **Coverage.** A machine check that every method the reference defines and can
  reach, and every CLI option of the three tools, has a counterpart in the port.
- **Real data.** Pfam-A.seed families, eight eggNOG5 Bacteria clusters and the
  truncated head of a ninth, and the mega-families of a Pfam-A.full window, each
  corpus pinned by the SHA-256 of its files.

### Scale validation

The four heaviest jobs — `strict`, `strictplus`, `automated1` and
`gappyout -htmlout` on a 3,583-sequence, ~7,400-column real EggNOG alignment
(O(sequence-count² × columns) conservation plus the HTML column walk) — are kept out
of the quick matrix and run on their own:

```bash
python3 verify/run_full.py --only scale
```

They are excluded from the quick gate purely for runtime; a timeout still counts as a
failure.

<!-- scale-table:start -->
| Job | Reference | Java | Speed-up |
|---|---:|---:|---:|
| `example.014.AA.EggNOG.COG0591.fasta -strict` | 57.1 s | 4.9 s | 12× |
| `example.014.AA.EggNOG.COG0591.fasta -strictplus` | 58.1 s | 5.1 s | 11× |
| `example.014.AA.EggNOG.COG0591.fasta -automated1` | 86.6 s | 7.7 s | 11× |
| `example.014.AA.EggNOG.COG0591.fasta -gappyout -htmlout` | 468.4 s | 6.3 s | 74× |
| `wide_60k.fa -sst` | 0.5 s | 0.1 s | 4× |
| `wide_60k.fa -sgt` | 0.5 s | 0.1 s | 5× |
| `wide_60k.fa -fasta` | 0.4 s | 0.1 s | 4× |
| `wide_60k.fa -phylip` | 0.6 s | 0.1 s | 7× |
| **total** | **672 s** | **24 s** | **28×** |
<!-- scale-table:end -->

Three changes deliver this while keeping every output byte identical: the
O(sequence-count² × columns) statistics loops run in parallel across rows/columns
with per-unit serial accumulation order (`-Dtrimal.threads=1` forces the plain serial
loop), the HTML writers precompute each column once per output block, and file-backed
writers batch their writes instead of issuing one flush syscall per line the way C++
`endl` does. The thread determinism gate is the proof that they change no byte.

### Mega-family scale (Pfam-A.full, up to 366,230 sequences)

The largest tier streams a 500 MB prefix window of **Pfam-A.full** (the 40.5 GB
official alignment archive) and extracts six mega-families — the largest, `7tm_1`
(GPCR), holds **366,230 sequences × 3,489 columns ≈ 1.28 billion cells**, 34× larger
than the biggest eggNOG cluster. Only O(n × L) methods apply here (the O(n²) identity
matrix is infeasible at this scale in *both* implementations):

<!-- mega-table:start -->
| Family — method | Reference | Java | Speed-up |
|---|---:|---:|---:|
| `pffull_01_7tm_1.fasta -gappyout` | 1513 s | 56 s | 27× |
| `pffull_01_7tm_1.fasta -noallgaps` | 1517 s | 60 s | 25× |
| `pffull_01_7tm_1.fasta -fasta` | 1520 s | 28 s | 54× |
| `pffull_01_7tm_1.fasta -sgc` | 1519 s | 30 s | 52× |
| `pffull_01_7tm_1.fasta -sgt` | 1542 s | 30 s | 51× |
| all 30 cases, up to 366K sequences | 189 min | 5 min | 38× |
<!-- mega-table:end -->

At this scale the reference spends 25 minutes and a sluggish per-column `std::string`
walk on `7tm_1` gappyout, while the port finishes in 56 seconds — and both write
identical 1.3 GB outputs.

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

Where the reference reads unallocated memory or dereferences a null pointer, the
port does something defined instead. These inputs are recorded as known divergences
with their exit status, streams and files kept, and are never counted as matches.

1. **Trimming a CLUSTAL, MEGA or PHYLIP input can crash the original.** Its `seqsInfo`
   pointer is NULL, but gap filtering, sequence selection, clustering and overlap
   filtering dereference it. The port uses empty auxiliary headers and goes on.
2. **Empty-alignment `-colnumbering` prints an arbitrary native integer.** The original
   reads beyond a zero-length array; the port prints an empty `#ColumnsMap` line. This
   narrowly defined difference is separately classified with rc/stderr/files checked.
3. **Malformed CRLF compare lists can trigger a native double-free.** The port retains
   the same failed-load message but returns a defined error rather than adding a
   glibc abort diagnostic. The companion suite records this separately.
4. **An empty comparison set, without `-forceselect`.** The original runs its
   comparison on a zero-element array and crashes; the port reports that the set
   lists no alignment and exits.
5. **A similarity matrix that did not load.** `trimal` ignores the failed load and
   dereferences the empty matrix as soon as a method needs it. The port reports the
   failure with the message `statal` uses for it and exits.
6. **`-htmlout` after a method returned nothing** (for example with contradictory
   manual boundaries). The original dereferences a null pointer; the port goes on as
   it does without `-htmlout`.
7. **Options with no input alignment, or with no sequence left to cluster.** The
   original indexes empty arrays: a consistency threshold, `-maxidentity`, a method
   with `-htmlout`. The port prints the reference's warning where the reference has
   one for the neighbouring case, and exits.
8. Zero-column slope, short-alignment and CDS edge cases avoid unsafe Java indexing
   where the source would perform undefined C++ memory accesses. The boundary audit
   records the actual native observations and the port's result without counting them
   as successful byte matches.
9. **Some of these crashes are not deterministic.** On six audit cases (a
   backtranslation with an empty selection) the reference segfaults in about half of
   the runs and exits normally in the rest. The port always takes the path the
   reference takes when it survives, and its output is byte-identical to it there.

There is no blanket normalization or SKIP rule hiding a reference-success versus
port-missing-output bug: a tagged vector is compared strictly on every run in which
the reference does not crash, and a difference on a vector without a tag fails its
gate.

## Project layout

```text
java/       Maven project, source, JUnit tests, recorded reference output and launchers
verify/     pinned reference metadata, fixtures, differential runners, fuzz gates and oracles
docs/       the method (method.md) and the initial source-archaeology reports
licenses/   retained third-party licenses
upstream/   local reference source/data snapshot; ignored by Git
```

## License

GPL-3.0-or-later; see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md). Adapted tracy-java
utilities retain their BSD-3-Clause notice. Ported glibc float algorithms retain
their LGPL-2.1-or-later/Sun notices. The corresponding licenses are included in the
JAR's `META-INF` directory. Distribution requires respecting the applicable source
and notice obligations; the reference executable is not a runtime dependency.
