#!/usr/bin/env python3
"""Generate deterministic small alignments and a valid-CLI differential matrix."""
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "verify/corpus/generated"
MANIFEST = ROOT / "verify/quick_manifest.txt"
SEED = 8151973


def fasta(name, records):
    path = CORPUS / name
    data = "".join(">" + identifier + "\n" + sequence + "\n"
                   for identifier, sequence in records)
    path.write_text(data, encoding="ascii", newline="\n")
    return path.relative_to(ROOT).as_posix()


def main():
    CORPUS.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    rows = []
    # Retain the existing real-data/error baseline, but keep the O(n^2 L) 3,583-row
    # conservation/strict jobs out of the quick suite. They remain in the full suite.
    baseline = ROOT / "verify/baseline_manifest.txt"
    excluded = []
    for line in baseline.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if "example.014." in line and any(
                flag in line for flag in (" -strict ", " -strictplus ", " -automated1 ", " -htmlout ")):
            excluded.append(line)
        else:
            rows.append(line)

    rows.extend(["--version", "-h"])

    reference = "upstream/dataset/example.001.AA.clw"
    small = "upstream/dataset/example.004.AA.fasta"
    real = "upstream/dataset/example.015.AA.bctoNOG.ENOG41099F3.fasta"
    rows.extend([
        f"-in {small} -nogaps -out $OUT",
        f"-in {small} -selectcols {{ 0-3,5,9-12 }} -out $OUT",
        f"-in {small} -selectseqs {{ 1-3 }} -out $OUT",
        f"-in {small} -sgt -out $OUT",
        f"-in {small} -sgc -sgt -ssc -sst -sident -soverlap -out $OUT",
        f"-in {small} -gt 0.8 -cons 60 -w 2 -colnumbering -out $OUT",
        f"-in {real} -ssc -sst -out $OUT",
        f"-in {real} -sgc -ssc -out $OUT",
        f"-in {small} -matrix upstream/dataset/matrix.BLOSUM62 -ssc -out $OUT",
        f"-compareset verify/corpus/cmp1.files -sfc -sft -out $OUT",
        f"-compareset verify/corpus/cmp1.files -cw 2 -sfc -out $OUT",
        f"-compareset verify/corpus/cmp1.files -w 2 -out $OUT",
        f"-compareset verify/corpus/cmp1.files -forceselect {reference} -sfc -out $OUT",
    ])

    # Fixed-seed AA and nucleotide inputs. No all-gap row/column is forced here; each
    # column has a real residue in the first row, avoiding undefined native paths.
    alphabets = ("ARNDCQEGHILKMFPSTWYV", "ACGT", "ACGTURYKMSWBDHV", "ACGU")
    fixtures = []
    for index in range(24):
        count = rng.choice([3, 5, 8, 12, 21])
        length = rng.choice([17, 31, 47, 60, 61, 101, 151, 219])
        alphabet = alphabets[index % len(alphabets)]
        template = "".join(rng.choice(alphabet) for _ in range(length))
        records = []
        for row in range(count):
            sequence = []
            for column, original in enumerate(template):
                if row == 0:
                    residue = original
                elif rng.random() < [0.08, 0.18, 0.32, 0.55][index % 4]:
                    residue = "-"
                elif rng.random() < 0.25:
                    residue = rng.choice(alphabet)
                else:
                    residue = original
                sequence.append(residue)
            records.append(("sample_%02d_%02d_long_name details:%d" % (index, row, row),
                            "".join(sequence)))
        fixture = fasta("random_%02d.fasta" % index, records)
        fixtures.append(fixture)
        for flags in ["-gappyout", "-strict", "-strictplus", "-automated1",
                      "-nogaps", "-noallgaps", "-gt 0.6", "-st 0.12",
                      "-gt 0.9 -cons 70", "-gt 0.8 -st 0.2 -cons 40 -w 2"]:
            rows.append(f"-in {fixture} {flags} -fasta -out $OUT")
        if index < 8:
            rows.extend([
                f"-in {fixture} -sgc -sgt -ssc -sst -out $OUT",
                f"-in {fixture} -sident -soverlap -out $OUT",
                ("SKIP:orig-ub-garbage-columnsmap " if index in (0, 4) else "")
                + f"-in {fixture} -gt 0.6 -complementary -colnumbering -out $OUT",
                f"-in {fixture} -gappyout -terminalonly -out $OUT",
                f"-in {fixture} -strict -block 3 -keepheader -out $OUT",
                f"-in {fixture} -gt 0.7 --set_boundaries {{ 2,{length - 3} }} -out $OUT",
                f"-in {fixture} -resoverlap 0.8 -seqoverlap 50 -out $OUT",
                f"-in {fixture} -clusters 3 -out $OUT",
                f"-in {fixture} -maxidentity 0.7 -out $OUT",
                f"-in {fixture} -gt 0.5 -htmlout $HTML -out $OUT",
                f"-in {fixture} -strict -htmlout $HTML -out $OUT",
            ])

    # Filename/header truncation and all output formats with nontrivial line wraps.
    for flags in ["-fasta", "-fasta_m10", "-clustal", "-nbrf", "-mega", "-nexus",
                  "-phylip", "-phylip_m10", "-phylip3.2", "-phylip3.2_m10",
                  "-phylip_paml", "-phylip_paml_m10"]:
        rows.append(f"-in {fixtures[7]} {flags} -out $OUT")
    rows.extend([
        f"-in {fixtures[1]} -fasta -keepheader -out $OUT",
        f"-in {fixtures[1]} -fasta_m10 -keepheader -out $OUT",
        f"-in {small} -out verify/corpus/generated/missing_parent/output.fasta",
    ])

    # Valid backtranslation inputs; ordinary, padding, stop-codon warning and errors.
    protein = fasta("protein.fasta", [("codingA", "MA-GFF"), ("codingB", "MALG-F")])
    cds = fasta("coding.fasta", [("codingA", "ATGGCTGGTTTTTTT"),
                                  ("codingB", "ATGGCTCTGGGTTTT")])
    stop = fasta("coding_stop.fasta", [("codingA", "ATGGCTTGATTTTTT"),
                                       ("codingB", "ATGGCTCTGGGTTTT")])
    rows.extend([
        f"-in {protein} -backtrans {cds} -out $OUT",
        f"-in {protein} -gt 0.5 -backtrans {cds} -out $OUT",
        f"-in {protein} -gt 0.5 -backtrans {cds} -mega -out $OUT",
        f"-in {protein} -backtrans {stop} -out $OUT",
        f"-in {protein} -backtrans {stop} -ignorestopcodon -out $OUT",
        f"-in {protein} -backtrans {stop} -splitbystopcodon -out $OUT",
    ])

    # Parameter error paths are separate from accepted arguments, including order.
    for arguments in ["-gt", "-st", "-ct", "-w", "-block", "-gt -1", "-st 1.1",
                      "-cons 101", "-nogaps -strict", "-strict -gt 0.3",
                      "-gt 0.3 -strict", "-w 2 -gappyout", "-gappyout -w 2",
                      "-st abc", "-resoverlap 0.2", "-seqoverlap 20",
                      "-matrix verify/corpus/missing.matrix", "-splitbystopcodon",
                      "-ignorestopcodon", "-colnumbering", "-complementary",
                      "--alternative_matrix unknown", "-keepheader -keepheader",
                      "-selectcols { 0,1000 }", "-selectseqs { 999 }",
                      "-gt 0.9 -sgc", "-fasta -clustal"]:
        rows.append(f"-in {small} {arguments}")

    MANIFEST.write_text("# Generated by verify/generate_suite.py, seed=%d\n" % SEED
                        + "\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    scale_rows = list(excluded) + generate_extended()
    (ROOT / "verify/scale_manifest.txt").write_text(
        "# Slow O(n^2 L) real-data jobs excluded from the quick suite.\n"
        + "\n".join(scale_rows) + "\n", encoding="utf-8", newline="\n")
    print("GENERATED cases=%d seed=%d" % (len(rows), SEED))
    print("QUICK " + str(MANIFEST))
    print("SCALE cases=%d" % len(scale_rows))
    generate_extended()


EXTENDED_SEED = 424242
EXT = ROOT / "verify/corpus/extended"
EXT_MANIFEST = ROOT / "verify/extended_manifest.txt"


def write_ext(name, text):
    EXT.mkdir(parents=True, exist_ok=True)
    path = EXT / name
    path.write_bytes(text.encode("latin1"))
    return path.relative_to(ROOT).as_posix()


def generate_extended():
    """Second, wider differential matrix. Independent seed; the quick manifest and its
    corpus stay byte-identical because this generator draws from its own Random."""
    rng = random.Random(EXTENDED_SEED)
    rows = []

    def add(line):
        rows.append(line)

    # ---------------- corpus ----------------
    aa = "ARNDCQEGHILKMFPSTWYV"
    dna = "ACGT"
    rna = "ACGU"
    deg = "ACGTURYKMSWBDHV"

    def make(alphabet, nseq, ncol, gap_rates, seed_offset, names=None, gapless=False):
        local = random.Random((EXTENDED_SEED, seed_offset).__hash__())
        records = []
        for r in range(nseq):
            seq = []
            for c in range(ncol):
                if gapless:
                    seq.append(alphabet[c % len(alphabet)] if r == 0 else
                               alphabet[(c * 7 + r * 3) % len(alphabet)])
                    continue
                if r == 0:
                    seq.append(alphabet[(c * 5) % len(alphabet)])
                elif local.random() < gap_rates[min(r, len(gap_rates) - 1)]:
                    seq.append("-")
                else:
                    seq.append(alphabet[local.randrange(len(alphabet))])
            name = names[r] if names else ("seq_%02d_of_long_name%d,tag:%d" % (r, r, r))
            records.append((name, "".join(seq)))
        return records

    def fasta_of(records):
        return "".join(">" + n + "\n" + s + "\n" for n, s in records)

    ext_aa = fasta_of(make(aa, 8, 60, [0.05, 0.12, 0.2, 0.3], 1))
    ext_dna = fasta_of(make(dna, 6, 47, [0.1, 0.25, 0.4], 2))
    ext_rna = fasta_of(make(rna, 5, 31, [0.15, 0.3], 3))
    ext_deg = fasta_of(make(deg, 7, 53, [0.08, 0.2, 0.35], 4))
    ext_aa2 = fasta_of(make(aa, 12, 77, [0.1, 0.22, 0.32], 5))
    edge_1seq = fasta_of(make(aa, 1, 30, [0.0], 6))
    edge_2x2 = fasta_of(make(dna, 2, 2, [0.0], 7))
    edge_1col = fasta_of(make(aa, 4, 1, [0.3], 8))
    edge_nogap = fasta_of(make(aa, 5, 40, [0.0], 9, gapless=True))
    edge_gaphalf = fasta_of(make(aa, 6, 50, [0.5, 0.6], 10))
    # identical gap counts everywhere: representative clustering hits its tie sort
    tie_seq = ["ACGT" * 11] * 9
    edge_ties = fasta_of([("tie_%02d" % i, s if i % 3 else s[:-1] + "-")
                          for i, s in enumerate(tie_seq)])
    edge_dup = fasta_of([("dup", "ACGTACGTAC"), ("other", "ACG--GTAC"),
                         ("dup", "ACGT---TAC"), ("last", "ACGTACG--")])
    edge_names = fasta_of([
        ("  leading_space_name", "ACGTACGTAC"),
        ("\ttab_name", "ACGT---TAC"),
        ("comma,name", "ACGTACG--"),
        ("colon:name", "ACGTACGTA"),
        ("a" * 80, "ACGTA-GTA"),
        ("ab;misleading", "ACGTT-GTA"),
    ])
    # format-detection quirks
    ext_pir = "".join(">P1;%s\nsome description %d\n%s*\n\n" % (n, i, s)
                      for i, (n, s) in enumerate(make(aa, 5, 40, [0.1, 0.2], 11)))
    deg_recs = make(deg, 5, 45, [0.1, 0.25], 12)
    ext_nex = ("#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX=5 NCHAR=45;\n"
               "FORMAT DATATYPE=DNA MISSING=? MATCHCHAR=. GAP=-;\n[matrix follows]\nMATRIX\n"
               + "".join(n + " " + s + "\n" for n, s in deg_recs) + ";\nEND;\n")
    mega_recs = make(dna, 4, 40, [0.1, 0.2], 13)
    ext_mega = ("#MEGA\n!Title extended title;\n!Format DataType=DNA;\n\n"
                + "".join("#" + n + "\n" + s + "\n\n" for n, s in mega_recs))
    ext_mega_seq = ("#mega\nTitle lowercase sequential;\n\n"
                    + "".join("#" + n + "\n" + s + "\n" for n, s in mega_recs))
    clw_recs = make(aa, 6, 48, [0.1, 0.2], 14)
    ext_clw = ("CLUSTAL multiple sequence alignment\n\n"
               + "".join(n + "  " + s + "\n" for n, s in clw_recs)
               + "\n      ** *   *  \n")
    ext_clw_lower = ext_clw.replace("CLUSTAL", "clustal", 1)
    ext_mislead = fasta_of([("ab;looks_pir", "ACGTACGTAC"), ("plain", "ACG---GTAC")])
    ext_crlf = fasta_of(make(aa, 4, 33, [0.15, 0.3], 15)).replace("\n", "\r\n")
    ext_nolf = fasta_of(make(aa, 4, 33, [0.15, 0.3], 16)).rstrip("\n")

    aa_path = write_ext("ext_aa.fa", ext_aa)
    dna_path = write_ext("ext_dna.fa", ext_dna)
    rna_path = write_ext("ext_rna.fa", ext_rna)
    deg_path = write_ext("ext_deg.fa", ext_deg)
    aa2_path = write_ext("ext_aa2.fa", ext_aa2)
    p1 = write_ext("edge_1seq.fa", edge_1seq)
    p2x2 = write_ext("edge_2x2.fa", edge_2x2)
    p1col = write_ext("edge_1col.fa", edge_1col)
    pnogap = write_ext("edge_nogap.fa", edge_nogap)
    phalf = write_ext("edge_gaphalf.fa", edge_gaphalf)
    pties = write_ext("edge_ties.fa", edge_ties)
    pdup = write_ext("edge_dupnames.fa", edge_dup)
    pnames = write_ext("edge_names.fa", edge_names)
    ppir = write_ext("ext_in.pir", ext_pir)
    pnex = write_ext("ext_in.nex", ext_nex)
    pmeg = write_ext("ext_in.meg", ext_mega)
    pmegseq = write_ext("ext_in_megaseq.meg", ext_mega_seq)
    pclw = write_ext("ext_in.clw", ext_clw)
    pclwlow = write_ext("ext_in_lower.clw", ext_clw_lower)
    pmis = write_ext("ext_mislead.fa", ext_mislead)
    pcrlf = write_ext("ext_crlf.fa", ext_crlf)
    pnolf = write_ext("ext_nolf.fa", ext_nolf)

    # custom matrices (symbols always cover the fixture alphabet)
    def matrix_text(symbols, labelled, asymmetric, seed):
        local = random.Random(seed)
        n = len(symbols)
        vals = [[local.choice([0, 0, 1, 2, -1, 4]) for _ in range(n)] for _ in range(n)]
        if not asymmetric:
            for i in range(n):
                for j in range(i + 1, n):
                    vals[i][j] = vals[j][i] = (vals[i][j] + vals[j][i]) // 2
        lines = [" \t ".join(symbols)]
        for i in range(n):
            lines.append(((symbols[i] + " ") if labelled else "")
                         + " ".join(str(v) for v in vals[i]))
        return "\n".join(lines) + "\n"

    m_aa = write_ext("matrix_aa.txt", matrix_text(aa, True, False, 21))
    m_asym = write_ext("matrix_asym_dna.txt", matrix_text("ACGT", False, True, 22))
    m_lower = write_ext("matrix_lower_dna.txt", matrix_text("acgt", True, False, 23))
    dna_fix = write_ext("matrix_dna.fa", fasta_of(make(dna, 5, 36, [0.15, 0.3], 24)))
    m_dna_fix = dna_fix

    # ---------------- full input x output conversion grid ----------------
    dataset_inputs = ["upstream/dataset/example.001.AA.clw",
                      "upstream/dataset/example.001.AA.msl",
                      "upstream/dataset/example.001.AA.phy",
                      "upstream/dataset/example.002.AA.phy",
                      "upstream/dataset/example.006.AA.pir",
                      "upstream/dataset/example.094.DNADeg.sequential_phy"]
    outputs = ["-clustal", "-fasta", "-fasta_m10", "-nbrf", "-nexus", "-mega",
               "-phylip", "-phylip_m10", "-phylip3.2", "-phylip3.2_m10",
               "-phylip_paml", "-phylip_paml_m10"]
    for input_file in dataset_inputs:
        for flag in outputs:
            add(f"-in {input_file} {flag} -out $OUT")
    for generated in [ppir, pnex, pmeg, pmegseq, pclw, pclwlow, pmis, pcrlf, pnolf,
                      pnames, pdup]:
        for flag in ["-fasta", "-clustal", "-phylip", "-nexus"]:
            add(f"-in {generated} {flag} -out $OUT")
    add(f"-in {ppir} -nbrf -out $OUT")            # PIR->PIR passthrough headers
    add(f"-in {pclw} -clustal -out $OUT")          # clustal title passthrough
    add(f"-in {aa_path} -keepheader -fasta -out $OUT")
    add(f"-in {aa_path} -keepheader -fasta_m10 -out $OUT")

    # ---------------- trimming method x parameter grids ----------------
    methods = ["-nogaps", "-noallgaps", "-gappyout", "-strict", "-strictplus", "-automated1"]
    # cleanStrict reads neighbourhood indices that are out of bounds below five
    # columns (upstream alignment.cpp:2665-2678): native survives by reading heap
    # garbage, which cannot be replicated deterministically, so the strict family is
    # only emitted for fixtures with at least five columns.
    strict_safe = [aa_path, dna_path, deg_path, pnogap, phalf, pties]
    trimmable = strict_safe + [p2x2, p1col]
    for fixture in trimmable:
        for method in (methods if fixture in strict_safe else
                       ["-nogaps", "-noallgaps", "-gappyout"]):
            add(f"-in {fixture} {method} -out $OUT")
    for fixture in [aa_path, dna_path]:
        for extra in ["-complementary", "-colnumbering", "-complementary -colnumbering",
                      "-keepseqs", "-terminalonly"]:
            for method in methods:
                add(f"-in {fixture} {method} {extra} -out $OUT")

    # manual threshold grids
    for value in ["0", "0.05", "0.2", "0.5", "0.8", "0.95", "1"]:
        add(f"-in {aa_path} -gt {value} -out $OUT")
    for value in ["0.05", "0.3", "0.7", "0.95"]:
        add(f"-in {aa_path} -st {value} -out $OUT")
        add(f"-in {aa_path} -gt 0.7 -st {value} -out $OUT")
    for value in ["0", "30", "60", "90", "100"]:
        add(f"-in {aa_path} -gt 0.6 -cons {value} -out $OUT")
        add(f"-in {aa_path} -st 0.2 -cons {value} -out $OUT")

    # block sizes across methods that accept them
    for size in ["1", "2", "4", "8", "20"]:
        for method in ["-gappyout", "-strict", "-gt 0.7"]:
            add(f"-in {aa_path} {method} -block {size} -out $OUT")

    # windows on methods and statistics
    for size in ["1", "2", "3", "7", "15"]:
        for method in ["-gt 0.7", "-gappyout", "-strict"]:
            add(f"-in {aa_path} {method} -w {size} -out $OUT")
        add(f"-in {aa_path} -sgc -w {size} -out $OUT")
        add(f"-in {aa_path} -sgt -gw {size} -out $OUT")
        add(f"-in {aa_path} -sst -sw {size} -out $OUT")
    add(f"-in {aa_path} -sgc -sgt -ssc -sst -out $OUT")
    add(f"-in {aa_path} -sst -ssc -sgt -sgc -out $OUT")
    add(f"-in {dna_path} -sgc -ssc -out $OUT")
    add(f"-in {rna_path} -sident -out $OUT")
    add(f"-in {deg_path} -soverlap -out $OUT")
    add(f"-in {aa_path} -w 2 -gw 3 -out $OUT")     # conflict error path
    add(f"-in {aa_path} -sw 2 -w 3 -out $OUT")

    # terminal boundaries
    add(f"-in {aa_path} -gappyout --set_boundaries {{ 0,59 }} -out $OUT")
    add(f"-in {aa_path} -gappyout --set_boundaries {{ 10,30 }} -out $OUT")
    add(f"-in {aa_path} -gt 0.7 --set_boundaries {{ 30,59 }} -terminalonly -out $OUT")
    add(f"-in {aa_path} -gt 0.7 --set_boundaries {{ 40,10 }} -out $OUT")
    add(f"-in {aa_path} -gt 0.7 --set_boundaries {{ 0,60 }} -out $OUT")

    # column / sequence selection grids
    for selection in ["0", "59", "0-59", "5,3,1", "0,2-4,10,58-59", "2-4"]:
        add(f"-in {aa_path} -selectcols {{ {selection} }} -out $OUT")
        add(f"-in {aa_path} -selectcols {{ {selection} }} -complementary -out $OUT")
    for selection in ["0-59", "60", "-1", "b"]:
        add(f"-in {aa_path} -selectcols {{ {selection} }} -out $OUT")
    for selection in ["0", "0-7", "1,3,5", "2-6"]:
        add(f"-in {aa2_path} -selectseqs {{ {selection} }} -out $OUT")
    add(f"-in {aa2_path} -selectseqs {{ 0,12 }} -out $OUT")
    add(f"-in {aa2_path} -selectseqs {{ 3-1 }} -out $OUT")

    # representative clustering incl. tie-heavy and full value range
    for k in range(1, 13):
        add(f"-in {aa2_path} -clusters {k} -out $OUT")
    for value in ["0.1", "0.3", "0.5", "0.7", "0.9"]:
        add(f"-in {aa2_path} -maxidentity {value} -out $OUT")
    add(f"-in {pties} -clusters 4 -out $OUT")
    add(f"-in {pties} -maxidentity 0.62 -out $OUT")
    # unaligned input + clustering segfaults the reference (null saveSequences write);
    # Java applies the documented safe-identity divergence instead.
    add("SKIP:orig-segv " + f"-in {pdup} -clusters 3 -out $OUT")

    # overlap filtering grid
    for res in ["0.3", "0.6", "0.9"]:
        for seq in ["10", "60", "100"]:
            add(f"-in {aa_path} -resoverlap {res} -seqoverlap {seq} -out $OUT")
    add(f"-in {phalf} -resoverlap 0.5 -seqoverlap 50 -keepseqs -out $OUT")

    # ---------------- trimmed output in every format ----------------
    for method in ["-gappyout", "-gt 0.7", "-strict"]:
        for flag in outputs:
            add(f"-in {aa_path} {method} {flag} -out $OUT")

    # ---------------- consistency comparisons ----------------
    cmp_a = write_ext("cmp_a.fa", ext_aa)
    reordered = fasta_of(list(reversed(make(aa, 8, 60, [0.05, 0.12, 0.2, 0.3], 1))))
    cmp_b = write_ext("cmp_b.fa", reordered)
    cmp_c = write_ext("cmp_c.fa", fasta_of(make(aa, 8, 52, [0.1, 0.25], 31)))
    cmp_d = write_ext("cmp_d.fa", fasta_of(make(aa, 8, 44, [0.2, 0.35], 32)))
    cmp_e = write_ext("cmp_e.fa", fasta_of(make(aa, 8, 66, [0.05, 0.18], 33)))
    cmp3 = write_ext("cmp3.files", "\n".join([cmp_a, cmp_b, cmp_c]) + "\n")
    cmp5 = write_ext("cmp5.files", "\n".join([cmp_a, cmp_b, cmp_c, cmp_d, cmp_e]) + "\n")
    cmp_first = write_ext("cmp_first.files", "\n".join([cmp_b, cmp_a, cmp_c]) + "\n")
    for list_file in [cmp3, cmp5, cmp_first]:
        add(f"-compareset {list_file} -out $OUT")
        add(f"-compareset {list_file} -sfc -sft -out $OUT")
        add(f"-compareset {list_file} -cw 2 -sfc -out $OUT")
        add(f"-compareset {list_file} -cw 3 -out $OUT")
        add(f"-compareset {list_file} -w 2 -out $OUT")
        for threshold in ["0.1", "0.3", "0.6", "0.9"]:
            add(f"-compareset {list_file} -ct {threshold} -out $OUT")
        add(f"-compareset {list_file} -ct 0.4 -cons 30 -out $OUT")
        add(f"-compareset {list_file} -ct 0.4 -complementary -out $OUT")
        add(f"-compareset {list_file} -ct 0.4 -htmlout $HTML -out $OUT")
    for forced in [cmp_a, cmp_b, cmp_c]:
        add(f"-compareset {cmp3} -forceselect {forced} -out $OUT")
        add(f"-compareset {cmp3} -forceselect {forced} -sfc -out $OUT")
        add(f"-compareset {cmp5} -forceselect {forced} -cw 2 -out $OUT")

    # ---------------- custom and alternative matrices ----------------
    add(f"-in {m_dna_fix} -matrix {m_asym} -st 0.3 -out $OUT")
    add(f"-in {m_dna_fix} -matrix {m_lower} -ssc -out $OUT")
    add(f"-in {m_dna_fix} -matrix {m_asym} -strict -out $OUT")
    add(f"-in {dna_path} --alternative_matrix degenerated_nt_identity -sgc -out $OUT")
    add(f"-in {deg_path} --alternative_matrix degenerated_nt_identity -sst -out $OUT")
    add(f"-in {dna_path} --alternative_matrix degenerated_nt_identity -automated1 -out $OUT")
    add(f"-in {m_dna_fix} -matrix {m_asym} -st 0.2 -htmlout $HTML -out $OUT")

    # ---------------- backtranslation grids ----------------
    protein = write_ext("bt_protein.fa", fasta_of([
        ("protB", "MA-GFFL"), ("protA", "MALGFF-"), ("protC", "M-AGFLI")]))
    cds_ok = write_ext("bt_cds_ok.fa", fasta_of([
        ("protA", "ATGGCTCTGGGCTTCTTT"), ("protB", "ATGGCTGGCTTCTTTCTG"),
        ("protC", "ATGGCTGGCTTTCTGAT")]))
    cds_taa = write_ext("bt_cds_taa.fa", fasta_of([
        ("protA", "ATGGCTTAAGGCTTCTTT"), ("protB", "ATGGCTGGCTTCTTTCTG"),
        ("protC", "ATGGCTGGCTTTCTGAT")]))
    cds_tag = write_ext("bt_cds_tag.fa", fasta_of([
        ("protA", "ATGGCTTAGGGCTTCTTT"), ("protB", "ATGGCTGGATTCTTTCTG"),
        ("protC", "ATGGCTGGCTTTCTGAT")]))
    cds_short = write_ext("bt_cds_short.fa", fasta_of([
        ("protA", "ATGGCT"), ("protB", "ATGGCTGGCTTCTTTCTG"),
        ("protC", "ATGGCTGGCTTTCTGAT")]))
    cds_long = write_ext("bt_cds_long.fa", fasta_of([
        ("protA", "ATGGCTCTGGGCTTCTTTAAA"), ("protB", "ATGGCTGGCTTCTTTCTG"),
        ("protC", "ATGGCTGGCTTTCTGAT")]))
    cds_gap = write_ext("bt_cds_gap.fa", fasta_of([
        ("protA", "ATGGCT--GGCTTC"), ("protB", "ATGGCTGGCTTCTTTCTG"),
        ("protC", "ATGGCTGGCTTTCTGAT")]))
    for cds in [cds_ok, cds_taa, cds_tag]:
        add(f"-in {protein} -backtrans {cds} -out $OUT")
        add(f"-in {protein} -backtrans {cds} -ignorestopcodon -out $OUT")
        add(f"-in {protein} -backtrans {cds} -splitbystopcodon -out $OUT")
        add(f"-in {protein} -gt 0.5 -backtrans {cds} -out $OUT")
        add(f"-in {protein} -gappyout -backtrans {cds} -nbrf -out $OUT")
        add(f"-in {protein} -gt 0.8 -backtrans {cds} -mega -out $OUT")
    add(f"-in {protein} -backtrans {cds_short} -out $OUT")
    add(f"-in {protein} -backtrans {cds_long} -out $OUT")
    add(f"-in {protein} -backtrans {cds_gap} -out $OUT")
    add(f"-in {protein} -gt 0.5 -backtrans {cds_gap} -out $OUT")

    # ---------------- miscellaneous option orders ----------------
    add(f"-in {aa_path} -out $OUT")
    add(f"-in {aa_path} -sgc -gappyout -out $OUT")
    add(f"-in {aa_path} -gappyout -sgc -out $OUT")
    add(f"-in {aa_path} -automated1 -sgc -ssc -colnumbering -out $OUT")
    add(f"-in {aa_path} -noallgaps -keepseqs -complementary -colnumbering -out $OUT")
    add(f"-in {phalf} -nogaps -keepseqs -out $OUT")
    add(f"-in {p1} -gappyout -out $OUT")
    add(f"-in {p1} -gt 0.5 -out $OUT")
    add(f"-in {p1} -htmlout $HTML -out $OUT")

    EXT_MANIFEST.write_text(
        "# Generated by verify/generate_suite.py, extended seed=%d\n" % EXTENDED_SEED
        + "\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print("EXTENDED cases=%d seed=%d" % (len(rows), EXTENDED_SEED))
    print("EXTENDED " + str(EXT_MANIFEST))

    # ---------------- ultra-wide alignment for the scale manifest ----------------
    # 4 x 60,000 columns: exercises the sorting-based statistics (quicksort over
    # 60k values) and deep loader paths. Deterministic; lives in the scale manifest
    # because native runtimes grow superlinearly on these.
    local = random.Random(EXTENDED_SEED + 1)
    wide_rows = []
    records = []
    for r in range(4):
        seq = "".join(local.choice("ACGT") if local.random() > 0.12 else "-"
                      for _ in range(60000))
        records.append(("wide_%d" % r, seq))
    wide = write_ext("wide_60k.fa", fasta_of(records))
    for flags in ["-sst", "-sgt", "-fasta", "-phylip"]:
        wide_rows.append(f"-in {wide} {flags} -out $OUT")
    print("WIDE cases=%d" % len(wide_rows))
    return wide_rows


if __name__ == "__main__":
    main()
