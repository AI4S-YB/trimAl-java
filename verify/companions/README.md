# readAl / statAl Java 11 companion verification

## Final result

The Java CLI implementations are direct transcriptions of `upstream/source/readAl.cpp`, `statAl.cpp`, and `defines.h`. They reuse the existing alignment, statistics, comparison, and `CppOut` implementations without resetting shared stdout formatting state.

| Suite | Cases | Byte-identical PASS | Actual FAIL | Known native UB divergence |
| --- | ---: | ---: | ---: | ---: |
| WSL Java 11 vs original native readal | 238 | 238 | 0 | 0 |
| WSL Java 11 vs original native statal | 237 | 236 | 0 | 1 |
| Combined WSL suite | 475 | 474 | 0 | 1 |
| Windows Java main entrypoints vs original Linux native | 49 | 49 | 0 | 0 |

The known native UB case is **not PASS**. Its `passed` field remains `false`; raw native/Java stderr, stdout, return codes, hashes, byte offsets, and artifact paths are preserved in `result.json`. The WSL script exits successfully only when there are no actual failures; `nonmatching` and `known_divergence` remain separate from the PASS count.

The final WSL run compiled all current Java sources with `javac --release 11 -encoding UTF-8` into `verify/companions/classes`, then packaged `verify/companions/companions.jar`. WSL runtime was OpenJDK 11.0.32.1. Both new CLI classes have class-file major version 55. Windows smoke used the installed OpenJDK 21.0.11 Java executable with the same Java-11-compatible jar; it is not claimed to be a Windows Java 11 runtime test.

Final independently compiled jar SHA-256:

`f0f29f08b7e1fab5719c3ff31b0e976856b8e79559fb3892177dbedbd7abac06`

Native tools report `v1.5.rev1 build[2025-11-25]`.

## Implemented entrypoints

- `java -cp <jar> trimal.ReadAl ...`
- `java -cp <jar> trimal.StatAl ...`

Both entrypoints restore the C `argv[0]` slot before invoking the transcribed parser. C main return value `-1` is mapped to process status `255` explicitly on Windows as well as Linux. Help, version, error strings, whitespace, upstream spelling mistakes, and the UTF-8 author bytes passed through the Latin-1 writer are preserved.

`ReadAl` supports all upstream output formats and name-shortening variants, format/type/info reports, reverse, unaligned `-onlyseqs`, HTML, stdout/file routing, immediate input loading, duplicate-option rejection, and source-order postchecks. It also translates Java file-open exceptions into the original output-file error.

`StatAl` supports all eight statistics flags (`-sgc`, `-sgt`, `-ssc`, `-sst`, `-sident`, `-scolidentt`, `-sfc`, `-sft`), `-i`/`-in`, comparison-set reference selection, forced reference, general/specific windows, custom and alternative matrices, source-order incompatibility checks, and the original output headings/order. Compare-list parsing emulates `std::getline` by retaining CR and not adding a synthetic EOF line.

## Coverage and artifacts

The 475-case WSL suite checks stdout, stderr, return code, and generated file bytes independently, without normalization. It includes:

- 30 deterministic independent random small datasets, seed 151, plus DNA/RNA/degenerate nucleotide/protein/aligned/unaligned/single/short/indeterminate/gap fixtures.
- Every readal format to stdout and file, every reverse-format combination, short names, only sequences, HTML, combined reports with output files, and native-generated input-format roundtrips.
- Help/no-arguments/version and their aliases or deliberate unsupported aliases; missing/repeated/unknown options; absent/empty/malformed inputs; output-directory and missing-parent failures; mutual-exclusion and postcheck errors.
- Every statistic alone, all six single-alignment statistics together and in reversed argument order, all eight with comparison sets, forced-reference and consistency windows.
- General/gap/similarity/consistency windows, numeric parsing quirks, too-large windows, malformed/custom/alternative matrices, and comparison-set loading/datatype/alignment errors.
- CRLF alignments, no-final-LF inputs and lists, raw UTF-8 and Latin-1 sequence names, case preservation, question-mark/dot punctuation, all-gap rows/columns, and repo-relative paths.
- Windows file byte comparisons, successful operations, and error exits explicitly observed as `255`.

Outputs are stored under `verify/companions/outputs/` and `verify/companions/windows_outputs/`. Final structured reports are `result.json` and `result_windows.json`; each report links cases to their raw artifacts and preserves byte counts and SHA-256 values.

## Known native UB difference

Case `statal:bytes_compareset_crlf` supplies a CRLF comparison list. Upstream `std::getline` retains CR, so the path ending in `dna.fa\r` fails to load. It then frees uninitialized comparison-list entries during cleanup (`upstream/source/statAl.cpp`, lines 305-310, 365-366, 468 in zero-based Read output).

Observed final native behavior is SIGABRT (Python return code `-6`) with its ordinary 92-byte alignment-load error followed by 32 bytes:

`double free or corruption (out)\n`

Java emits the identical 92-byte error prefix and exits `255`, with empty stdout and no output file on both sides. The script classifies this only when the exact test name, native signal, Java failure status, identical stdout/files, identical error prefix, CR-bearing filename, and exact glibc diagnostic all match. Any timeout, changed byte output, Java success, or different failure remains an actual FAIL. This native crash is not fabricated in the portable Java implementation.

## Initial failure and shared-layer resolution

The first suite produced 417 cases: 416 PASS and one actual FAIL. An amino-acid input with `--alternative_matrix degenerated_nt_identity` reached `SimilarityMatrix.java:257` and accessed the 15x15 identity table past its bounds after allocating a 20x20 protein matrix. Native emitted 855 bytes with status 0; Java emitted an exception with status 1. The shared numeric-module owner reproduced the reference binary's contiguous-table layout and corrected that implementation. No CLI-side fallback or skip was added.

In the final run `statal:alternative_aa` is byte-identical PASS, stdout length 855 and SHA-256 `c43c98727e4acc095991301cbb71eb3b9f777c5840e350bd4f62f888cbe128c1` on both sides. `history_initial.log` preserves the initial 416/417 failure log. `history_bytes.log` and `result_bytes_.json` preserve the earlier 57/58 boundary run before strict known-UB classification was introduced; they are historical reports, not the final result.

## Re-running

From the repository root in WSL/Linux:

```sh
python3 verify/companions/verify_companions.py
```

This rebuilds **only** the companion-specific classes/jar and does not invoke Maven, the main harness, or `verify/work`.

Reuse the independent jar:

```sh
python3 verify/companions/verify_companions.py --no-compile
```

Test a coordinator-produced integrated jar without any compilation:

```sh
python3 verify/companions/verify_companions.py --java-jar java/target/trimal.jar
```

`--java-jar` accepts the actual relative or absolute integrated jar path and automatically skips compilation. `--filter <substring>` can run a focused case subset; it writes a separate `result_<filter>.json` rather than replacing the complete final report.

Windows smoke runs inside WSL and launches the Windows Java executable. Set the executable to the installed Java path when needed:

```sh
python3 verify/companions/verify_windows_companions.py --java-exe '/mnt/c/Program Files/Eclipse Adoptium/jdk-21.0.11.10-hotspot/bin/java.exe'
```

It also accepts `--java-jar <integrated-jar>` and performs no compilation. Roots are calculated from each script's location; CLI data paths and report artifact paths are repository-relative. All generated inputs, outputs, logs, classes, jars, and reports stay under `verify/companions`.

## Scope

Only `java/src/main/java/trimal/ReadAl.java`, `StatAl.java`, and companion-specific verification files were created or edited by this task. Existing alignment/statistics/utilities, Maven output, the main harness, `verify/work`, git state, and remote repositories were not modified.
