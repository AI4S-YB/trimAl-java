#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Deterministic CLI argument fuzz gate for trimAl.

Builds seed-fixed argument vectors out of the real option names and pools of values
(in range, out of range, non-numeric, empty, huge, duplicated, contradictory, with
the value missing at the end) and runs the strict differential harness over them.
Byte equality of rc/stderr/stdout/files is required for every vector, whether the
reference accepts it or rejects it, so the argument parser and its error messages
are compared as strictly as the trimming paths. Run on Linux/WSL:
    python3 verify/argfuzz_suite.py                 # generate + run (label "argfuzz")
    python3 verify/argfuzz_suite.py --all           # every label
    python3 verify/argfuzz_suite.py --generate-only

The inputs are small alignments (tens of KB at most), so even the O(n^2 x L)
methods of the reference finish quickly. Each label is a fixed seed with its own
manifest (verify/<label>_manifest.txt); the harness reports how many vectors the
reference accepted (`VACUITY ... reference_ok=`), so a label cannot degenerate into
an argument-error-only test unnoticed.
"""
import argparse
from pathlib import Path
import random
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

# label -> (seed, number of argument vectors). "argfuzz3" draws its vectors from
# VALID_KINDS: the shapes of commands the reference is known to accept, with random
# values and inputs. The random labels are rejected by the argument validation about
# five times out of six; this one reaches the trimming paths.
LABELS = {"argfuzz": (7001, 400), "argfuzz2": (7002, 400), "argfuzz3": (7003, 400)}
VALID_LABELS = {"argfuzz3"}

# Vectors (1-based positions in the manifest) on which the reference build crashes
# with a signal. On the random labels: selection lists out of range, an empty
# comparison set, a similarity matrix that did not load, a method that returned
# nothing while -htmlout is set, ... On the accepted-shape label (argfuzz3) it is
# ordinary trimming that crashes: a method on a CLUSTAL, MEGA or PHYLIP input whose
# seqsInfo is NULL, an empty result, a selection or overlap filter on such an input.
# They carry the harness's SKIP:orig-segv tag, which counts them as KNOWN_DIVERGENCE
# only while the reference really crashes; should it stop crashing on one of them,
# that vector turns into an ordinary byte comparison.
NATIVE_CRASH = {
    "argfuzz": {2, 7, 30, 47, 48, 50, 76, 92, 102, 121, 161, 169, 208, 215, 232, 245,
                246, 287, 292, 310, 325, 327, 347, 367, 371, 381},
    "argfuzz2": {14, 18, 34, 58, 60, 61, 84, 89, 112, 126, 139, 140, 144, 151, 171,
                 208, 230, 273, 291, 292, 307, 320, 323, 379},
    "argfuzz3": {1, 7, 18, 36, 50, 54, 63, 66, 74, 75, 79, 83, 96, 109, 112, 125, 143,
                 164, 181, 190, 209, 212, 226, 239, 241, 243, 248, 257, 269, 277, 301,
                 302, 311, 337, 350, 365, 376, 383, 387},
}

INPUTS = [
    "upstream/dataset/example.001.AA.clw", "upstream/dataset/example.001.AA.msl",
    "upstream/dataset/example.001.AA.phy", "upstream/dataset/example.002.AA.clw",
    "upstream/dataset/example.002.AA.phy", "upstream/dataset/example.003.AA.clw",
    "upstream/dataset/example.004.AA.fasta", "upstream/dataset/example.005.AA.fasta",
    "upstream/dataset/example.006.AA.pir", "upstream/dataset/example.007.AA.fasta",
    "upstream/dataset/example.009.AA.fasta", "upstream/dataset/example.010.AA.fasta",
    "upstream/dataset/example.092.DNA.fasta", "upstream/dataset/example.093.DNA.fasta",
    "upstream/dataset/example.094.DNADeg.sequential_phy",
    "upstream/dataset/example.055.AA.bctoNOG.ENOG4109GY9.fasta",
    "upstream/dataset/example.043.AA.bctoNOG.ENOG4109ADN.fasta",
    "upstream/dataset/alignments_comparison.1", "upstream/dataset/alignments_comparison.2",
    "upstream/dataset/alignments_comparison.3",
    "verify/corpus/extended/ext_aa.fa", "verify/corpus/extended/ext_dna.fa",
    "verify/corpus/extended/ext_in.pir", "verify/corpus/extended/bt_protein.fa",
    "verify/corpus/one_seq.fasta", "verify/corpus/bad_trunc.fasta",
    "verify/corpus/empty.fasta",
]
# an input path that deliberately does not exist
ABSENT = "verify/corpus/missing.fasta"
METHODS = ["-gappyout", "-strict", "-strictplus", "-automated1", "-nogaps", "-noallgaps"]
STATS = ["-sgc", "-sgt", "-ssc", "-sst", "-sident", "-soverlap", "-sfc", "-sft"]
FORMATS = ["-fasta", "-fasta_m10", "-clustal", "-phylip", "-phylip_m10", "-phylip3.2",
           "-phylip3.2_m10", "-phylip_paml", "-phylip_paml_m10", "-nexus", "-nbrf", "-mega"]
MODIFIERS = ["-complementary", "-colnumbering", "-terminalonly", "-keepseqs", "-keepheader"]
FLOAT_OPTIONS = ["-gt", "-gapthreshold", "-st", "-simthreshold", "-ct", "-conthreshold",
                 "-maxidentity", "-resoverlap"]
INT_OPTIONS = ["-w", "-gw", "-sw", "-cw", "-block", "-seqoverlap"]
FLOATS = ["0", "0.3", "0.5", "0.7", "0.9", "1", "1.0", "1.0000001", "-0.1", "-1", "2",
          "100", "1e-3", ".5", "abc", "", "nan", "inf", "0.3333333333333333",
          "0.9999999999999999", "1e400", " 0.5", "0.", "-.5"]
INTS = ["0", "1", "2", "3", "5", "10", "-1", "-5", "3.7", "abc", "", "2147483647",
        "2147483648", "99999999999", "0x10", "1e2", "+2", " 4"]
NUMBERS = ["0", "25", "50", "60", "99", "100", "101", "-1", "1.5", "abc", ""]
SELECTIONS = ["0", "0-3", "5,3,1", "0,0,0", "5-1", "999", "-1", "b", "", "0-3,5,9-12", "1, 2", "0,", ","]
BOUNDARIES = ["0,59", "10,30", "40,10", "0", "0,1,2", "abc", "-1,5", "", "5,5", "1,2147483647"]
MATRICES = ["upstream/dataset/matrix.BLOSUM62", "upstream/dataset/matrix.Degenerated_DNA",
            "verify/corpus/missing.matrix", "upstream/dataset/example.004.AA.fasta"]
ALTERNATIVES = ["0", "1", "2", "15", "16", "-1", "abc", "1000000", ""]
UNKNOWN = ["-phylip3.2", "-nogap", "-automated2", "-colscore", "-identity", "-wrongflag",
           "-1", "-overlap", "-mw", "-IN", "-", "--"]
BACKTRANS = ["verify/corpus/extended/bt_cds_ok.fa", "verify/corpus/extended/bt_cds_gap.fa",
             "verify/corpus/extended/bt_cds_short.fa", "verify/corpus/extended/bt_cds_long.fa",
             "verify/corpus/extended/bt_cds_taa.fa", "verify/corpus/extended/bt_cds_tag.fa",
             "verify/corpus/generated/coding.fasta", "verify/corpus/generated/coding_stop.fasta"]


def pick(rng, pool, low, high):
    """A random subset of `pool` with between low and high members, order kept."""
    return rng.sample(pool, rng.randrange(low, high + 1))


def groups(rng, pool, low, high):
    """The same for a pool of option groups (each a list of tokens), flattened."""
    return [token for group in pick(rng, pool, low, high) for token in group]


def out(rng, probability=0.85):
    tokens = []
    if rng.random() < probability:
        tokens += ["-out", "$OUT"]
    if rng.random() < 0.12:
        tokens += ["-htmlout", "$HTML"]
    return tokens


def value_option(rng):
    option = rng.choice(FLOAT_OPTIONS + INT_OPTIONS + ["-cons"])
    pool = FLOATS if option in FLOAT_OPTIONS else INTS if option in INT_OPTIONS else NUMBERS
    return [option, rng.choice(pool)]


def k_plain(rng):
    return ["-in", rng.choice(INPUTS)] + pick(rng, METHODS + FORMATS + MODIFIERS, 1, 3) + out(rng)


def k_threshold(rng):
    tokens = ["-in", rng.choice(INPUTS)]
    for _ in range(rng.randrange(1, 4)):
        tokens += value_option(rng)
    return tokens + pick(rng, METHODS + STATS + MODIFIERS, 0, 2) + out(rng)


def k_window(rng):
    tokens = ["-in", rng.choice(INPUTS)]
    for option in pick(rng, ["-w", "-gw", "-sw", "-cw"], 1, 4):
        tokens += [option, rng.choice(INTS)]
    return tokens + pick(rng, STATS, 1, 3) + out(rng)


def k_select(rng):
    tokens = ["-in", rng.choice(INPUTS), rng.choice(["-selectcols", "-selectseqs"]), "{",
              rng.choice(SELECTIONS), "}"]
    if rng.random() < 0.15:  # no braced block: the reference rejects the argument
        tokens = tokens[:3] + [tokens[4]]
    return tokens + pick(rng, MODIFIERS + FORMATS, 0, 2) + out(rng)


def k_boundaries(rng):
    method = rng.choice([["-gappyout"], ["-strict"], ["-gt", rng.choice(FLOATS[:8])], []])
    return (["-in", rng.choice(INPUTS)] + method + ["--set_boundaries", "{",
            rng.choice(BOUNDARIES), "}"] + pick(rng, ["-terminalonly", "-colnumbering"], 0, 1)
            + out(rng))


def k_cluster(rng):
    return (["-in", rng.choice(INPUTS), "-clusters", rng.choice(INTS)]
            + groups(rng, [["-maxidentity", rng.choice(FLOATS[:10])], ["-keepseqs"],
                           ["-colnumbering"]], 0, 2) + out(rng))


def k_overlap(rng):
    return (["-in", rng.choice(INPUTS), "-resoverlap", rng.choice(FLOATS),
             "-seqoverlap", rng.choice(INTS)] + pick(rng, ["-keepseqs", "-keepheader"], 0, 1)
            + out(rng))


def compare_lists(rng, count, directory):
    """Small lists of input paths: valid, with a missing file, blank lines or CRLF.
    They live in a directory of the label's own: the manifests refer to them by path,
    so one label's lists must never overwrite another's."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for index in range(count):
        members = pick(rng, INPUTS[:20] + [ABSENT], 1, 3)
        text = "\n".join(members) + "\n"
        style = rng.randrange(5)
        if style == 1:
            text = text.replace("\n", "\r\n")
        elif style == 2:
            text = "\n" + text + "\n"
        elif style == 3:
            text = text.rstrip("\n")
        elif style == 4:
            text = ""
        path = directory / ("set_%02d.files" % index)
        path.write_bytes(text.encode("ascii"))
        paths.append(path.relative_to(ROOT).as_posix())
    return paths


def k_compare(rng, lists):
    tokens = ["-compareset", rng.choice(lists)] + pick(rng, ["-sfc", "-sft"], 0, 2)
    for _ in range(rng.randrange(0, 3)):
        tokens += value_option(rng)
    if rng.random() < 0.3:
        tokens += ["-forceselect", rng.choice(INPUTS[:20])]
    return tokens + out(rng)


def k_backtrans(rng):
    extra = groups(rng, [["-ignorestopcodon"], ["-splitbystopcodon"],
                         ["-gt", rng.choice(FLOATS[:8])], ["-mega"], ["-fasta_m10"]], 0, 2)
    return (["-in", rng.choice(["verify/corpus/extended/bt_protein.fa",
                                "verify/corpus/generated/protein.fasta"] + INPUTS[:8]),
             "-backtrans", rng.choice(BACKTRANS)] + extra + out(rng))


def k_matrix(rng):
    tokens = ["-in", rng.choice(INPUTS)]
    if rng.random() < 0.7:
        tokens += ["-matrix", rng.choice(MATRICES)]
    else:
        tokens += ["--alternative_matrix", rng.choice(ALTERNATIVES)]
    return tokens + pick(rng, STATS, 1, 3) + out(rng)


def k_stats(rng):
    return ["-in", rng.choice(INPUTS)] + pick(rng, STATS, 1, 4) + out(rng, 0.5)


def k_noinput(rng):
    return pick(rng, METHODS + STATS + MODIFIERS, 1, 3) + (value_option(rng) if rng.random() < 0.5 else [])


def k_duplicate(rng):
    option = rng.choice(["-in", "-out", "-gt", "-w", "-fasta"])
    pool = {"-in": INPUTS[:8], "-out": ["$OUT", "$HTML", "verify/corpus/argfuzz"],
            "-gt": FLOATS[:8], "-w": INTS[:8], "-fasta": ["-fasta"]}[option]
    tokens = ["-in", rng.choice(INPUTS), option, rng.choice(pool), option, rng.choice(pool)]
    return tokens + out(rng, 0.5)


def k_conflict(rng):
    return (["-in", rng.choice(INPUTS)] + pick(rng, METHODS, 2, 2) + pick(rng, FORMATS, 0, 2)
            + out(rng, 0.5))


def k_dangling(rng):
    option = rng.choice(FLOAT_OPTIONS + INT_OPTIONS + ["-in", "-out", "-htmlout", "-matrix",
                                                       "-compareset", "-backtrans", "-selectcols"])
    tokens = ["-in", rng.choice(INPUTS)] if rng.random() < 0.6 and option != "-in" else []
    return tokens + pick(rng, METHODS + MODIFIERS, 0, 2) + [option]


def k_help(rng):
    return rng.choice([["-h"], ["--version"], ["-h", "-in", rng.choice(INPUTS)],
                       ["--version", "-h"], ["-h", "-h"], ["--version", "extra"],
                       ["-in", rng.choice(INPUTS), "-h"], ["-v"], ["--help"]])


def k_unknown(rng):
    tokens = ["-in", rng.choice(INPUTS)] if rng.random() < 0.6 else []
    return tokens + [rng.choice(UNKNOWN)] + pick(rng, METHODS, 0, 1)


KINDS = [
    (k_plain, 14), (k_threshold, 12), (k_window, 8), (k_select, 8), (k_boundaries, 5),
    (k_cluster, 5), (k_overlap, 4), (k_compare, 10), (k_backtrans, 6), (k_matrix, 6),
    (k_stats, 5), (k_noinput, 4), (k_duplicate, 4), (k_conflict, 4), (k_dangling, 5),
    (k_help, 3), (k_unknown, 3),
]


# ---- shapes of commands the reference accepts, with random values and inputs ----

VALID_FLOATS = ["0", "0.1", "0.3", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]
VALID_INTS = ["1", "2", "3", "4", "5", "8", "10", "20"]
VALID_CONS = ["0", "25", "40", "50", "60", "75", "99", "100"]
VALID_SELECTS = ["0", "0-3", "5,3,1", "0-3,5,9-12", "0-20", "10-5"]
VALID_BOUNDARIES = ["0,59", "10,30", "0,1,2", "5,5", "2,20"]
# the real alignment files among INPUTS (not the empty, single-sequence or broken ones)
ALIGNMENTS = INPUTS[:20]
PROTEINS = [path for path in ALIGNMENTS if ".AA." in path]
CODING = ["verify/corpus/extended/bt_cds_ok.fa", "verify/corpus/extended/bt_cds_gap.fa",
          "verify/corpus/extended/bt_cds_taa.fa", "verify/corpus/generated/coding.fasta",
          "verify/corpus/generated/coding_stop.fasta"]
# one alignment in three formats: every member loads, and the reference compares them
COMPARABLE = [
    ["upstream/dataset/example.001.AA.clw", "upstream/dataset/example.001.AA.msl",
     "upstream/dataset/example.001.AA.phy"],
]


def valid_lists(rng, count, directory):
    """Comparison lists whose members all exist and are the same alignment."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for index in range(count):
        members = sorted(rng.sample(COMPARABLE[index % len(COMPARABLE)], rng.randrange(2, 4)))
        path = directory / ("v_%02d.files" % index)
        path.write_text("\n".join(members) + "\n", encoding="ascii", newline="\n")
        paths.append(path.relative_to(ROOT).as_posix())
    return paths


def v_value(rng, option):
    pool = (VALID_CONS if option == "-cons" else
            VALID_FLOATS if option in FLOAT_OPTIONS else VALID_INTS)
    return [option, rng.choice(pool)]


def v_method(rng, lists):
    return (["-in", rng.choice(ALIGNMENTS), rng.choice(METHODS)]
            + pick(rng, FORMATS[:9] + ["-fasta_m10"] + MODIFIERS, 0, 2) + out(rng))


def v_threshold(rng, lists):
    tokens = ["-in", rng.choice(ALIGNMENTS)]
    for option in pick(rng, ["-gt", "-st", "-cons", "-w", "-gw", "-sw", "-block"], 1, 3):
        tokens += v_value(rng, option)
    return (tokens + pick(rng, ["-complementary", "-colnumbering", "-keepseqs", "-keepheader",
                                "-terminalonly"] + FORMATS[:9], 0, 2) + out(rng))


def v_window_stats(rng, lists):
    tokens = ["-in", rng.choice(ALIGNMENTS)]
    for option in pick(rng, ["-w", "-gw", "-sw", "-cw"], 0, 2):
        tokens += v_value(rng, option)
    return tokens + pick(rng, STATS, 1, 4) + out(rng, 0.4)


def v_select(rng, lists):
    return (["-in", rng.choice(ALIGNMENTS), rng.choice(["-selectcols", "-selectseqs"]), "{",
             rng.choice(VALID_SELECTS), "}"] + pick(rng, MODIFIERS + FORMATS[:9], 0, 2) + out(rng))


def v_boundaries(rng, lists):
    method = rng.choice([["-gappyout"], ["-strict"], ["-automated1"],
                         ["-gt", rng.choice(VALID_FLOATS)],
                         ["-selectcols", "{", rng.choice(VALID_SELECTS), "}"]])
    return (["-in", rng.choice(ALIGNMENTS)] + method + ["--set_boundaries", "{",
            rng.choice(VALID_BOUNDARIES), "}"] + (["-terminalonly"] if rng.random() < 0.7 else [])
            + out(rng))


def v_cluster(rng, lists):
    tokens = ["-in", rng.choice(ALIGNMENTS), "-clusters", rng.choice(VALID_INTS)]
    if rng.random() < 0.5:
        tokens += ["-maxidentity", rng.choice(VALID_FLOATS)]
    return tokens + pick(rng, ["-keepseqs", "-colnumbering"], 0, 1) + out(rng)


def v_overlap(rng, lists):
    return (["-in", rng.choice(ALIGNMENTS), "-resoverlap", rng.choice(VALID_FLOATS),
             "-seqoverlap", rng.choice(["10", "25", "40", "50", "60", "75"])]
            + pick(rng, ["-keepseqs", "-keepheader"], 0, 1) + out(rng))


def v_backtrans(rng, lists):
    protein = rng.choice(["verify/corpus/extended/bt_protein.fa",
                          "verify/corpus/generated/protein.fasta"])
    extra = groups(rng, [["-ignorestopcodon"], ["-splitbystopcodon"],
                         ["-gt", rng.choice(VALID_FLOATS)], ["-gappyout"]], 0, 1)
    return (["-in", protein, "-backtrans", rng.choice(CODING)] + extra
            + pick(rng, ["-mega", "-fasta", "-clustal", "-nexus"], 0, 1) + out(rng))


def v_matrix(rng, lists):
    return (["-in", rng.choice(PROTEINS), "-matrix", "upstream/dataset/matrix.BLOSUM62"]
            + pick(rng, ["-ssc", "-sst", "-sident", "-strict", "-automated1", "-strictplus"], 1, 2)
            + out(rng))


def v_compare(rng, lists):
    tokens = ["-compareset", rng.choice(lists)] + pick(rng, ["-sfc", "-sft"], 1, 2)
    for option in pick(rng, ["-w", "-cw"], 0, 1):
        tokens += v_value(rng, option)
    return tokens + out(rng)


VALID_KINDS = [
    (v_method, 12), (v_threshold, 14), (v_window_stats, 8), (v_select, 8), (v_boundaries, 8),
    (v_cluster, 7), (v_overlap, 6), (v_backtrans, 9), (v_matrix, 8), (v_compare, 10),
]


def quote(token):
    # braces and the harness placeholders are single tokens; quote only what shlex needs
    return token if token in ("{", "}", "$OUT", "$HTML") else shlex.quote(token)


def generate(label, seed, count):
    rng = random.Random(seed)
    absent = [path for path in INPUTS if not (ROOT / path).is_file()]
    if absent:
        raise SystemExit("missing inputs (run verify/generate_suite.py first): " + ", ".join(absent))
    valid = label in VALID_LABELS
    directory = ROOT / "verify/corpus" / label
    lists = valid_lists(rng, 8, directory) if valid else compare_lists(rng, 24, directory)
    kinds = VALID_KINDS if valid else KINDS
    rows = []
    crashes = NATIVE_CRASH.get(label, ())
    for position in range(1, count + 1):
        kind = rng.choices([kind for kind, _ in kinds], weights=[weight for _, weight in kinds])[0]
        tokens = kind(rng, lists) if valid or kind is k_compare else kind(rng)
        row = " ".join(quote(token) for token in tokens)
        rows.append(("SKIP:orig-segv " if position in crashes else "") + row)
    manifest = ROOT / ("verify/%s_manifest.txt" % label)
    manifest.write_text("# Generated by verify/argfuzz_suite.py, seed=%d\n" % seed
                        + "\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print("ARGFUZZ label=%s cases=%d seed=%d" % (label, len(rows), seed))
    return manifest


def run(manifest, timeout, jobs, jar):
    command = [sys.executable, str(ROOT / "verify/differential.py"), str(manifest),
               "--timeout", str(timeout), "--jobs", str(jobs)]
    if jar:
        command += ["--jar", jar]
    return subprocess.run(command).returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--label", choices=sorted(LABELS), default="argfuzz")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--jar", help="JAR to test (default: java/target/trimal.jar)")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--jobs", type=int, default=4)
    options = parser.parse_args()
    status = 0
    for label in (sorted(LABELS) if options.all else [options.label]):
        seed, count = LABELS[label]
        manifest = generate(label, seed, count)
        if not options.generate_only:
            print("RUN", label, flush=True)
            status = max(status, run(manifest, options.timeout, options.jobs, options.jar))
    return status


if __name__ == "__main__":
    sys.exit(main())
