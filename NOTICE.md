# Third-party provenance

The Java implementation is a derivative of **trimAl v1.5.1**, originally developed
by Victor Fernández-Rodríguez, Salvador Capella-Gutiérrez and Toni Gabaldón
(2009–2020). The upstream source identifies itself as
`trimAl v1.5.rev1 build[2025-11-25]`.

- Upstream: https://github.com/inab/trimal
- Primary license: GNU GPL version 3; see [LICENSE](LICENSE). Upstream source headers
  describe the last available GPL version. This rewrite is distributed under
  GPL-3.0-or-later, subject to third-party notices below.
- The reference build uses `g++ -Wall -O2` under Linux x86-64. The precise acquired
  source files are fingerprinted in `verify/reference-source-sha256.json`.

## Reused Java utilities

`CppFormat` (C++ iostream floating-point formatting) was initially adapted from the
locally available **tracy-java** project and has since been extended for trimAl's
variable precision and fixed-point output. Its original BSD-3-Clause notice,
including the EMBL and tracy-java authors' copyright statements, is retained in
[licenses/tracy-java-BSD-3-Clause.txt](licenses/tracy-java-BSD-3-Clause.txt).

Other components once carried over from tracy-java (a C++ writer, a Latin-1
helper, boost-style rounding and a libstdc++ introsort replica) were later removed
as trimAl's own transcribed utilities superseded them. trimAl's sorting algorithm
is its own last-element-pivot quicksort, transcribed iteratively with an identical
permutation.

## Single-precision math compatibility

`Utils.expf` and `Utils.log10f` include ports of glibc 2.35 float arithmetic code and
coefficient tables. The originating Free Software Foundation code (2017–2022)
uses LGPL-2.1-or-later; the full license is retained in
[licenses/LGPL-2.1.txt](licenses/LGPL-2.1.txt). The logarithm implementation also
retains the Sun Microsystems 1993 permissive notice in the source. Exact source
paths and original notices are embedded beside the ported algorithms. This code
is Java-only; it does not call a native library at runtime.

No original reference binary, generated test output, build artifact or downloaded
upstream archive is intended to be committed as part of the Java source project.
