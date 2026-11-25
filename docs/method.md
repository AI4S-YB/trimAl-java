# Method: behaviour first, driven by an oracle

How the port was built, and on what grounds its outputs are known to equal the
reference's. This is the record of the approach for anyone who has to reproduce it on
another program.

## 1. Pin the oracle before writing a line

The oracle is the reference executable itself, not its documentation and not the
algorithm as described in a paper:

- upstream tag `v1.5.1`, commit `d637091abe33595775f40480970d1a18d87a7bcb`;
- built with `g++ 11.4.0 -Wall -O2` on Ubuntu 22.04 x86-64 (glibc 2.35);
- all 24 source, header and makefile blobs verified against the tag's tree by Git
  blob hash (`verify/identify_reference.py`), in CI on every run.

The build environment is part of the oracle. The float overloads of `expf` and
`log10f` return bits that depend on the libm they come from, and the reference uses
them directly in its conservation statistic. The same C++ oracle that matches the
port on glibc 2.35 disagrees on 1,906 of 8,219 probes on glibc 2.39. So the
environment is recorded with every result, and CI runs on the same image.

"Equal" means all of: exit status, stdout, stderr, and every file written, compared
byte for byte. No whitespace, precision or line-ending normalisation is applied
anywhere. A single differing byte is a failure.

## 2. Transcribe, do not redesign

The port keeps the reference's shape so that reasoning about one is reasoning about
the other: one Java class per C++ source file, the same method names and
signatures, the same loop order, the same float and double promotion points, the
same unstable quicksort, the same output-stream state that persists between
functions.

Java has no counterpart for several things C++ code relies on. Each was found by
reading the source and by differential failures, and each is emulated by one small
unit with its own bit-level oracle:

| Reference behaviour | Emulation |
|---|---|
| `std::cout` formatting: `%g`, precision, `fixed`, `setw`, left/right flags that stick | `util/CppFormat`, `util/CppOut` |
| `std::string::operator[]` at `size()` yields `'\0'` | `Utils.charAt` |
| `strcpy` stops at the first NUL of a line read with `getline` | `Alignment.Cursor.readLine` |
| float `expf` / `log10f` from glibc (tables and polynomials) | `Utils.expf`, `Utils.log10f` |
| `atoi`, `atof`, `is_number` of the C library | `Utils.atoi`, `Utils.atof`, `Utils.isNumber` |
| conversion of NaN and out-of-range doubles to `int` on x86 | `Utils.x86IntCast` |
| `main` returning `-1`, which the shell reports as 255 | entry classes normalise the status |

Nothing is inferred from "what the tool should do". Where the reference prints a
misspelled message, a misleading number or a result that looks wrong, the port prints
the same bytes.

## 3. Prove equality from the outside

A differential harness runs the reference and the port on the same arguments in the
same directory and compares everything they produce (`verify/differential.py`). The
cases are listed in manifests; `$OUT`, `$HTML` and `$IN` expand to per-case paths.

Three rules keep the record honest:

1. **A difference is never hidden.** A case whose reference crashes, or prints an
   integer read past the end of an array, is recorded as *known divergence*. It is
   counted apart from the matches, its exit status, streams and files are kept, and
   the classification applies only while the reference really does crash: a tagged
   vector on which the reference exits normally is compared strictly.
2. **A suite must prove that the reference ran.** Every report states how many cases
   the reference exited 0 on. A suite on which it never does compares error messages
   only, and fails unless that is declared on the command line. This guard was added
   late and immediately showed that three gates were empty: a corpus that had been
   written without its gaps (the reference rejected every case as unaligned), and two
   manifests using options the companion tools do not have.
3. **Evidence names its build.** Results are only accepted from a run in a clean tree,
   and the record carries the commit, the JAR's SHA-256, the reference binaries'
   SHA-256 and versions, the toolchain, and the machine load during each gate.

## 4. Layers of evidence

Each layer can fail alone, and none is derived from another:

| Layer | What it pins down | Size |
|---|---|---|
| Module oracles | float arithmetic, matrices, conservation, gaps, formatting, parsing, quicksort, the two `pow()` sites | 8,219 + 492,300 bit-level probes |
| Differential matrices | the CLI on generated and reference inputs, every option, error paths | 1,415 cases |
| Loader fuzz | mutated real files: truncation, byte flips, CRLF, control bytes, rewrapping | 3 × 180 |
| Argument fuzz | random and accepted-shaped argument vectors | 3 × 400 |
| Boundary audit | edge shapes, every writer read back by every parser, direct method probes | 605 + 3 |
| Golden CLI gate | the port run in-process against recorded reference output, in `mvn test` | 1,432 cases, any OS and JDK |
| Real data | Pfam-A.seed families, eggNOG clusters (and the truncated head of one of them), Pfam-A.full mega-families, companions on them | 262 cases |
| Cross-platform | Windows Java against the Linux reference | 1,000+ cases |
| Determinism | 1, 2, 7 and default thread counts against each other and the reference | 4 workloads |
| Coverage | every method the reference can reach and every CLI option exists in the port | machine-checked |

The suites overlap; the counts are reported per suite and never summed into one
coverage figure.

## 5. Undefined behaviour in the reference

The reference performs operations that C++ leaves undefined. Each site was put in
one of three classes, and the class decides what the port does:

- **Deterministic in this build and observable.** The all-above-threshold
  conservation read, for example, yields a fixed value here. The port reproduces
  that value. This is an observation about one build, stated as such.
- **Garbage.** An integer read from unallocated memory changes between runs. The
  port prints a defined value instead, and the case is a known divergence: there is
  nothing to be equal to.
- **Crashes.** A null `seqsInfo` dereferenced when a CLUSTAL or PHYLIP input is
  trimmed, an empty comparison set, a similarity matrix that failed to load and is
  used anyway. The port either continues with a defined result or reports the
  problem and exits; it never prints a stack trace where the reference crashes.
  Some of these crashes are not deterministic: on six audit cases (a
  backtranslation with an empty selection) the reference segfaults in about half of
  the runs and exits normally in the rest. Whenever it exits normally, its output is
  byte-identical to the port's — the strongest evidence there is that the port's
  safe continuation is the right one. Such a case is a strict match in one run and a
  recorded native crash in the next, never a failure in either.

The class of every site is visible in the result reports. The port never turns an
undefined behaviour of the reference into a matching success.

## 6. Speed, only where bytes cannot change

The port is faster than the reference on large inputs because of three changes that
were each proven byte-neutral rather than argued to be so:

- the O(sequences² × columns) statistics run in parallel, with each thread owning
  its own output slots and accumulating in the same serial order as the loop it
  replaces;
- the HTML writers compute each column once per output block instead of once per
  row;
- file output is batched instead of flushing after every line.

The thread determinism gate runs the heavy workloads with 1, 2, 7 and the default
number of threads and requires identical bytes across all of them and against the
reference. `-Dtrimal.threads=1` selects the plain serial loop.

## 7. The record

`verify/run_full.py` runs every gate on one build and writes one record per gate;
`verify/validation_summary.py` turns the record into `verify/validation-summary.json`
and renders the README tables from it. A number in those tables is never typed by
hand. The summary refuses a run made in a dirty tree, a run that mixes commits or
JARs, and any run with a failure.
