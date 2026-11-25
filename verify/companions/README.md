# readAl / statAl companion verification

The `readal` and `statal` command lines (`trimal.ReadAl`, `trimal.StatAl`) are
transcriptions of `upstream/source/readAl.cpp`, `statAl.cpp` and `defines.h`. They
reuse the alignment, statistics, comparison and `CppOut` classes of the port. Both
entry points restore the C `argv[0]` slot before running the transcribed parser, and
map a `main` return value of `-1` to the process status `255` on every platform.

## Run

From the repository root in WSL or Linux:

```bash
python3 verify/companions/verify_companions.py --java-jar java/target/trimal.jar
python3 verify/companions/verify_windows_companions.py --java-jar java/target/trimal.jar
```

Both scripts test `java/target/trimal.jar` by default and skip compilation when
`--java-jar` names a JAR. `verify_companions.py` is the gate the `run_full.py` driver
calls `companions`; `verify_windows_companions.py` needs a Windows `java.exe`
reachable from WSL (`TRIMAL_WINDOWS_JAVA`, or on PATH) and is the driver's
`windows_companions` gate. `--filter <substring>` runs a subset of the Linux suite.

Each run writes its report to a new directory under `verify/companions/runs/`
(`result.json`, `result_windows.json`); the generated inputs, the outputs of every
case and the comparison lists stay under `verify/companions/` and are ignored by Git.
Every report links each case to its raw artifacts and keeps byte counts and SHA-256
values.

## What the Linux suite compares

475 cases, each on exit status, stdout, stderr and the produced file, with no
normalization:

- 30 deterministic random small datasets (seed 151), and DNA, RNA, degenerate
  nucleotide, protein, aligned, unaligned, single-sequence, short, indeterminate and
  gap fixtures.
- Every `readal` output format to stdout and to a file, every reverse/format
  combination, short names, `-onlyseqs`, HTML, combined reports, and input files
  that the native `readal` wrote itself.
- Help, no arguments, version and their aliases; missing, repeated and unknown
  options; absent, empty and malformed inputs; an output directory and a missing
  parent directory; mutual exclusions and the late checks.
- Every `statal` statistic alone, all six single-alignment statistics together and in
  reversed order, all eight with comparison sets, forced reference and consistency
  windows.
- General, gap, similarity and consistency windows, numeric parsing quirks, too-large
  windows, malformed, custom and alternative matrices, and comparison-set loading,
  datatype and alignment errors.
- CRLF alignments, inputs and lists without a final line feed, UTF-8 and Latin-1
  sequence names, case preservation, `?`/`.` punctuation, all-gap rows and columns.

## What the Windows suite compares

49 cases: the Windows JVM runs the same `trimal.ReadAl` / `trimal.StatAl` entry
points on the same inputs, and everything it prints and writes is compared with what
the Linux reference produced on the same arguments.

## Classification

A case is `pass` only when every byte matches. `known_native_ub` is a case in which
the reference itself fails undefined behaviour; it keeps `passed: false`, its raw
streams and hashes, and is never counted as a pass. The scripts exit non-zero only
on a `fail`, but report the three counts apart.

There is one such case. `statal:bytes_compareset_crlf` supplies a CRLF comparison
list. `std::getline` keeps the CR, so the path ending in `dna.fa\r` does not load,
and `statAl` then frees uninitialized comparison-list entries during cleanup
(`statAl.cpp`). The native process dies with SIGABRT after printing its ordinary
92-byte load error followed by `double free or corruption (out)`. Java prints the
identical 92 bytes and exits `255`, with empty stdout and no output file on either
side. The script classifies this only when the test name, the native signal, the Java
status, the streams, the files, the CR-bearing file name and the exact glibc
diagnostic all match; any timeout, changed output, Java success or other failure is
an actual `fail`.

## A bug this suite found

Its first version had one real failure: `statal` on an amino-acid input with
`--alternative_matrix degenerated_nt_identity` read the 15×15 identity table past
its bounds after allocating a 20×20 protein matrix (`SimilarityMatrix`). The native
program printed 855 bytes and exited 0; Java threw. The port now reproduces the
reference binary's contiguous-table layout, and the case is a byte-identical pass.
No fallback or skip was added.
