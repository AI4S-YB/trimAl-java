#!/usr/bin/env python3
"""Independent small byte-differential audit. Run with WSL Python; no Maven/shared work."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / 'verify/alignment_audit'
CORPUS = AUDIT / 'corpus'
SEED = 1511011


def write(name, content):
    path = CORPUS / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content.encode('latin1') if isinstance(content, str) else content)
    return str(path.relative_to(ROOT))


def formats(names, seqs):
    n, length = len(names), len(seqs[0])
    fasta = ''.join('>' + name + ' original header,tag:x\n' + seq + '\n' for name, seq in zip(names, seqs))
    clustal = 'CLUSTAL audit\n\n' + ''.join(name + '  ' + seq + '\n' for name, seq in zip(names, seqs))
    clustal_blocks = 'CLUSTAL audit\n\n'
    for start in range(0, length, 7):
        clustal_blocks += ''.join(name + '  ' + seq[start:start+7] + '  7\n' for name, seq in zip(names, seqs)) + '     *:*\n\n'
    pir = ''.join('>DL;' + name + '\n' + name + ' description\n' + seq + '*\n\n' for name, seq in zip(names, seqs))
    nexus = '#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX=%d NCHAR=%d;\nFORMAT DATATYPE=DNA MISSING=? MATCHCHAR=. GAP=-;\nMATRIX\n' % (n, length)
    nexus += ''.join(name + ' ' + seq + '\n' for name, seq in zip(names, seqs)) + ';\nEND;\n'
    nexus_blocks = '#nexus\nbegin data;\ndimensions ntax=%d nchar=%d;\nformat datatype=dna gap=-;\nmatrix\n[comment]\n' % (n, length)
    for start in range(0, length, 7):
        nexus_blocks += ''.join(name + ' ' + seq[start:start+7] + '\n' for name, seq in zip(names, seqs)) + '\n'
    nexus_blocks += ';\nend;\n'
    mega = '#MEGA\n!Title audit title;\n!Format DataType=DNA;\n\n' + ''.join('#' + name + '\n' + seq + '\n\n' for name, seq in zip(names, seqs))
    mega_blocks = '#mega\nTitle lowercase tag;\n\n'
    for start in range(0, length, 7):
        mega_blocks += ''.join('#' + name + ' ' + seq[start:start+7] + '\n' for name, seq in zip(names, seqs)) + '\n'
    phylip = '%d %d\n' % (n, length) + ''.join(name + ' ' + seq + '\n' for name, seq in zip(names, seqs))
    phylip_blocks = '%d %d\n' % (n, length)
    phylip_blocks += ''.join(name + ' ' + seq[:7] + '\n' for name, seq in zip(names, seqs)) + '\n'
    for start in range(7, length, 7):
        phylip_blocks += ''.join(seq[start:start+7] + '\n' for seq in seqs) + '\n'
    return dict(fasta=fasta, clustal=clustal, clustal_blocks=clustal_blocks, pir=pir,
                nexus=nexus, nexus_blocks=nexus_blocks, mega=mega, mega_blocks=mega_blocks,
                phylip=phylip, phylip_blocks=phylip_blocks)


def generate():
    rng = random.Random(SEED)
    cases = []
    def case(name, file, opts=(), group='parse', **extra):
        cases.append(dict(name=name, file=file, opts=list(opts), group=group, **extra))
    names = ['alpha01', 'beta02', 'gamma03', 'delta04', 'gaps05']
    seqs = ['ACGT--ACGTACNNACGTACGT--', 'ACGT--AC-TACNNACGTACGT--',
            'AC-T--ACGTACNNAC-TACGT--', 'ACGT--ACGTACNNACGTACGT--', '-' * 24]
    docs = formats(names, seqs)
    for fmt, text in docs.items():
        for variant, data in [('lf', text), ('crlf', text.replace('\n', '\r\n')),
                              ('no_final_lf', text.rstrip('\n')),
                              ('leading_spaces', '\n'.join('   ' + line if line else line for line in text.split('\n'))),
                              ('leading_tabs', '\n'.join('\t\t' + line if line else line for line in text.split('\n'))),
                              ('tabs_then_space', '\n'.join('\t ' + line if line else line for line in text.split('\n')))]:
            file = write(fmt + '_' + variant + '.txt', data)
            case(fmt + '_' + variant, file, ['-fasta'])
            if variant in ('lf', 'crlf', 'no_final_lf'):
                case(fmt + '_' + variant + '_strict', file, ['-strict', '-fasta', '-colnumbering'], 'trim')
        file = write(fmt + '_lf.txt', text)
        for out in ['-clustal', '-nbrf', '-fasta', '-phylip', '-phylip3.2', '-phylip_paml', '-nexus', '-mega', '-onlyseqs']:
            case(fmt + '_to_' + out[1:], file, [out], 'convert')
        for op in [['-nogaps'], ['-noallgaps'], ['-selectcols', '{', '0,2,8-10', '}'],
                   ['-selectseqs', '{', '0,4', '}'], ['-clusters', '2']]:
            case(fmt + '_op_' + op[0][1:], file, op + ['-fasta', '-keepheader', '-colnumbering'], 'trim')
    fasta = write('fasta_lf.txt', docs['fasta'])
    operations = [
        ['-nogaps'], ['-noallgaps'], ['-gappyout'], ['-strict'], ['-strictplus'], ['-automated1'],
        ['-gt', '0.9'], ['-gt', '0.99', '-cons', '80'], ['-gt', '0.99', '-cons', '100'],
        ['-st', '0.5', '-cons', '80'], ['-gt', '0.99', '-st', '0.5', '-cons', '80'],
        ['-selectcols', '{', '0,2,8-10', '}'], ['-selectseqs', '{', '0,4', '}'],
        ['-resoverlap', '0.8', '-seqoverlap', '75'], ['-clusters', '1'], ['-clusters', '2'],
        ['-clusters', '5'], ['-maxidentity', '0.8'], ['-nogaps', '-terminalonly'],
        ['-gt', '0.8', '-block', '3']]
    for idx, op in enumerate(operations):
        for flags in [[], ['-keepseqs'], ['-complementary'], ['-keepseqs', '-complementary']]:
            case('fasta_trim_%02d_%s' % (idx, '_'.join(x[1:] for x in flags) or 'normal'),
                 fasta, op + flags + ['-fasta', '-keepheader', '-colnumbering'], 'trim')
    for header in ['>a', '>ab', '>abc', '> a', '>\ta', '>a:x', '>a,b', '>', '> ', '>\t']:
        file = write('short_header_%02d.fa' % len(cases), header + '\nACGT\n>b\nACGT\n')
        case('short_header_%02d' % len(cases), file, ['-fasta', '-keepheader'], 'header')
    for name, data in [('empty', b''), ('blank', b' \t\n\r\n'), ('only_header', b'>a'),
                       ('nul', b'>alpha\nAC\x00GT\n>beta\nAC\x00GT\n'),
                       ('ascii_digits', b'>alpha\nAC1T\n>beta\nAC1T\n'),
                       ('high_byte', b'>alpha\nAC\xe9T\n>beta\nAC\xe9T\n'),
                       ('punctuation', b'>alpha\nAC._?!T\n>beta\nAC._?!T\n')]:
        file = write(name + '.fa', data)
        case('invalid_' + name, file, ['-fasta'], 'invalid')
    for n, length in [(1, 1), (2, 1), (2, 2), (3, 3), (5, 4), (7, 5), (9, 17), (31, 47), (23, 199)]:
        rows = [''.join(rng.choice('ACGTNN---') for _ in range(length)) for _ in range(n)]
        names2 = ['random%02d' % i for i in range(n)]
        file = write('seed_%d_%d.fa' % (n, length), formats(names2, rows)['fasta'])
        for idx, op in enumerate([['-strict'], ['-strictplus'], ['-gt', '0.95', '-cons', '80'],
                                 ['-st', '0.5', '-cons', '60'], ['-clusters', str(min(2, n))],
                                 ['-sident', '-soverlap']]):
            case('seed_%d_%d_%d' % (n, length, idx), file, op + ['-fasta'], 'seed')
    # GUI-independent I/O error tests use exactly the same literal file name.
    case('missing_input', 'verify/alignment_audit/corpus/does_not_exist.fa', ['-fasta'], 'io')
    case('output_missing_parent', fasta, ['-fasta', '-out', '{MISSING_OUT}'], 'io')
    case('html_missing_parent', fasta, ['-nogaps', '-fasta', '-htmlout', '{MISSING_HTML}'], 'io')
    case('output_directory', fasta, ['-fasta', '-out', 'verify/alignment_audit/corpus'], 'io')
    case('save_output', fasta, ['-nogaps', '-fasta', '-out', '{OUT}', '-keepheader', '-colnumbering'], 'io')
    case('save_html', fasta, ['-nogaps', '-fasta', '-out', '{OUT}', '-htmlout', '{HTML}', '-colnumbering'], 'io')
    # Backtranslation includes reordered CDS records, rare residues, and exact-end NUL.
    protein = write('protein.fa', '>alpha\nM-AQK\n>beta\nMA-QK\n')
    cds = write('cds.fa', '>beta\nATGGCTCAGAAA\n>alpha\nATGGCTCAGAAA\n')
    mega_cds = write('cds_title.meg', '#MEGA\n!Title retained CDS title;\n\n#beta\nATGGCTCAGAAA\n\n#alpha\nATGGCTCAGAAA\n\n')
    for idx, op in enumerate([[], ['-nogaps'], ['-selectcols', '{', '0,2', '}'],
                             ['-selectseqs', '{', '0', '}'], ['-clusters', '1']]):
        for output in ['-fasta', '-mega', '-nexus']:
            case('backtrans_%d_%s' % (idx, output[1:]), protein,
                 op + ['-backtrans', mega_cds, output], 'backtrans')
    for idx, (prot, dna) in enumerate([('MAQK', 'ATGTGACAGAAA'), ('MUQK', 'ATGTGACAGAAA'),
                                      ('mOxb', 'ATGTGATAATAG'), ('MAQK', 'ATGTAACAGAAA'),
                                      ('MAQK', 'ATGTAGCAGAAA'), ('MAQK', 'ATGtgaCAGAAA'),
                                      ('MAQK', 'ATGGCTCAGAA'), ('MAQK', 'ATGGCTCAGAAATGA'),
                                      ('MAQK', 'ATGTGATAATAG'), ('MAQK', 'ATGGCTCAGAAATGATGA'),
                                      ('MAQK', 'ATG-GTCAGAAA'), ('MAXX', 'ATGGCT')]):
        pp = write('stop_protein_%d.fa' % idx, '>alpha\n' + prot + '\n>beta\n' + prot + '\n')
        dd = write('stop_cds_%d.fa' % idx, '>alpha\n' + dna + '\n>beta\n' + dna + '\n')
        for flag in [[], ['-ignorestopcodon'], ['-splitbystopcodon']]:
            case('stop_%d_%s' % (idx, flag[0][1:] if flag else 'normal'), pp,
                 ['-backtrans', dd, '-fasta'] + flag, 'stop')
    # Defined-domain trimming: no sequence made wholly of gaps, common retained cores.
    core_names = ['alpha01', 'beta02', 'gamma03', 'delta04']
    core_rows = ['--ACGTACGTACGTACGTACGT--', 'A-ACGTACGTACGTACGTACGT--',
                 '--ACGTACGT-CGTACGTACGTAA', '--ACGTACGTACGTACGT-CGT--']
    valid_ops = [['-nogaps'], ['-noallgaps'], ['-strict'], ['-strictplus'], ['-gappyout'],
                 ['-gt', '0.99', '-cons', '80'], ['-st', '0.5', '-cons', '80'],
                 ['-gt', '0.99', '-st', '0.5', '-cons', '80'], ['-clusters', '2'],
                 ['-selectcols', '{', '0,2,8-10', '}'], ['-selectseqs', '{', '0', '}'],
                 ['-resoverlap', '0.8', '-seqoverlap', '75'], ['-nogaps', '-terminalonly']]
    for fmt, text in formats(core_names, core_rows).items():
        file = write('core_' + fmt + '.txt', text)
        for idx, op in enumerate(valid_ops):
            case('core_%s_trim_%02d' % (fmt, idx), file,
                 op + ['-fasta', '-keepheader', '-colnumbering'], 'defined_trim')
        if fmt == 'fasta':
            for output in ['-nbrf', '-mega', '-nexus', '-clustal', '-phylip', '-phylip3.2']:
                case('core_roundtrip_' + output[1:], file, [output], 'roundtrip', roundtrip=True)
    # C++ strict intentionally uses seqsInfo[0..newN), not the selected sequence rows.
    prefix_bug = write('strict_unfiltered_info.fa', formats(['removed00'] + core_names,
                                                          ['-' * 24] + core_rows)['fasta'])
    for op in ['-strict', '-strictplus']:
        case('strict_unfiltered_info_' + op[1:], prefix_bug,
             [op, '-fasta', '-keepheader', '-colnumbering'], 'defined_trim')
    for length in [10, 50, 60, 101]:
        file = write('writer_boundary_%d.fa' % length,
                     formats(['alpha01', 'beta02'], ['ACGT' * (length // 4) + 'ACGT'[:length % 4]] * 2)['fasta'])
        for output in ['-nbrf', '-mega', '-nexus']:
            case('writer_boundary_%d_%s' % (length, output[1:]), file, [output], 'roundtrip', roundtrip=True)
    core_fasta = write('core_fasta.txt', formats(core_names, core_rows)['fasta'])
    case('core_save_output', core_fasta, ['-nogaps', '-fasta', '-out', '{OUT}', '-keepheader', '-colnumbering'], 'io')
    case('core_save_html', core_fasta, ['-nogaps', '-fasta', '-out', '{OUT}', '-htmlout', '{HTML}', '-colnumbering'], 'io')
    case('backtrans_empty_selected_cols', protein, ['-selectcols', '{', '0-4', '}', '-backtrans', cds, '-fasta'], 'backtrans')
    case('backtrans_empty_selected_seqs', protein, ['-selectseqs', '{', '0-1', '}', '-backtrans', cds, '-fasta'], 'backtrans')
    case('backtrans_keepseqs_empty_cols', protein, ['-selectcols', '{', '0-4', '}', '-keepseqs', '-backtrans', cds, '-fasta'], 'backtrans')
    for count in [1, 2, 3, 5]:
        rows = ['A' * count, 'A' * count]
        file = write('small_protein_%d.fa' % count, formats(['alpha01', 'beta02'], rows)['fasta'])
        dnafile = write('small_cds_%d.fa' % count, formats(['alpha01', 'beta02'], ['GCT' * count] * 2)['fasta'])
        for op in [['-strict'], ['-strictplus'], ['-nogaps']]:
            case('small_backtrans_%d_%s' % (count, op[0][1:]), file, op + ['-backtrans', dnafile, '-fasta'], 'backtrans')
    (AUDIT / 'cases.json').write_text(json.dumps(dict(seed=SEED, cases=cases), indent=2), encoding='utf8')
    return cases


def run_one(case, classes, results):
    name = case['name']
    result_dir = results / name
    result_dir.mkdir(parents=True, exist_ok=True)
    relative = str(result_dir.relative_to(ROOT))
    replacements = {'{OUT}': relative + '/shared.out', '{HTML}': relative + '/shared.html',
                    '{MISSING_OUT}': relative + '/missing/child.out', '{MISSING_HTML}': relative + '/missing/child.html'}
    opts = [replacements.get(arg, arg) for arg in case['opts']]
    args = ['-in', case['file']] + opts
    if case.get('roundtrip'):
        args += ['-out', replacements['{OUT}']]
    record = dict(case, args=args)
    artifacts = [ROOT / replacements['{OUT}'], ROOT / replacements['{HTML}']]
    values = []
    for impl, prefix in [('cpp', ['/usr/bin/timeout', '--signal=KILL', '3s', str(ROOT / 'verify/bin/trimal')]),
                         ('java', ['/usr/bin/java', '-XX:-UsePerfData', '-cp', str(classes.relative_to(ROOT)), 'trimal.TrimAl'])]:
        for artifact in artifacts:
            if artifact.exists():
                artifact.unlink()
        try:
            proc = subprocess.run(prefix + args, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=6)
            rc, out, err = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as e:
            rc, out, err = 124, e.stdout or b'', e.stderr or b''
        (result_dir / (impl + '.stdout')).write_bytes(out)
        (result_dir / (impl + '.stderr')).write_bytes(err)
        files = []
        for idx, artifact in enumerate(artifacts):
            value = artifact.read_bytes() if artifact.is_file() else None
            files.append(value)
            if value is not None:
                (result_dir / (impl + ('.out' if idx == 0 else '.html'))).write_bytes(value)
        stages = []
        if case.get('roundtrip') and rc == 0 and files[0] is not None:
            reread_args = ['-in', replacements['{OUT}'], '-fasta']
            try:
                reread = subprocess.run(prefix + reread_args, cwd=ROOT, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, timeout=6)
                stages = [reread.returncode, reread.stdout, reread.stderr]
            except subprocess.TimeoutExpired as e:
                stages = [124, e.stdout or b'', e.stderr or b'']
            (result_dir / (impl + '.reread.stdout')).write_bytes(stages[1])
            (result_dir / (impl + '.reread.stderr')).write_bytes(stages[2])
        values.append((rc, out, err, files, stages))
        record[impl] = dict(rc=rc, stdout_sha=hashlib.sha256(out).hexdigest(), stderr_sha=hashlib.sha256(err).hexdigest(),
                            exception=b'Exception' in err or (bool(stages) and b'Exception' in stages[2]),
                            output_present=[v is not None for v in files], reread_rc=stages[0] if stages else None)
    cpp, java = values
    record['match'] = cpp == java
    record['diff'] = [label for idx, label in enumerate(['rc', 'stdout', 'stderr', 'files', 'reread']) if cpp[idx] != java[idx]]
    record['status'] = 'match' if record['match'] else 'mismatch'
    record['roundtrip_executed'] = bool(case.get('roundtrip') and all(v[0] == 0 and v[4] and v[4][0] == 0 for v in values))
    # Classification never changes the raw comparison or turns a divergence into pass.
    if record['match']:
        record['classification'] = 'exact'
    elif cpp[0] in (-9, 124, 137) or java[0] in (-9, 124, 137):
        record['classification'] = 'timeout'
    elif cpp[0] in (-11, -6, 139):
        record['classification'] = 'native_crash_safety_divergence'
    elif (record['diff'] == ['stdout'] and cpp[0] == 0 and cpp[2] == java[2]
          and cpp[3] == java[3] and java[1] == b'#ColumnsMap\t\n'
          and cpp[1].startswith(b'#ColumnsMap\t')
          and cpp[1][len(b'#ColumnsMap\t'):].strip().lstrip(b'-').isdigit()):
        record['classification'] = 'known_native_ub_empty_map'
        record['source_evidence'] = 'alignment.cpp:1807-1809 reads saveResidues[0] when residNumber=0'
    else:
        record['classification'] = 'actual_fail'
    (result_dir / 'result.json').write_text(json.dumps(record, indent=2), encoding='utf8')
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--round', default='round1')
    parser.add_argument('--filter', default='')
    parser.add_argument('--generate-only', action='store_true')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.environ['LC_ALL'] = 'C'
    CORPUS.mkdir(parents=True, exist_ok=True)
    cases = generate()
    if args.generate_only:
        print('generated', len(cases), 'cases, seed', SEED)
        return
    classes = AUDIT / ('classes_' + args.round)
    results = AUDIT / ('results_' + args.round)
    classes.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    sources = sorted((ROOT / 'java/src/main/java').rglob('*.java'))
    proc = subprocess.run(['/usr/bin/javac', '--release', '11', '-encoding', 'UTF-8', '-d', str(classes)] + [str(p) for p in sources],
                          cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
    (results / 'javac.stdout').write_bytes(proc.stdout)
    (results / 'javac.stderr').write_bytes(proc.stderr)
    if proc.returncode:
        sys.stdout.buffer.write(proc.stderr)
        raise SystemExit(proc.returncode)
    (results / 'source_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}, indent=2))
    cases = [c for c in cases if args.filter in c['name']]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(lambda c: run_one(c, classes, results), cases))
    summary = dict(seed=SEED, count=len(records), match=sum(r['match'] for r in records),
                   cpp_success=sum(r['cpp']['rc'] == 0 for r in records),
                   cpp_success_match=sum(r['cpp']['rc'] == 0 and r['match'] for r in records),
                   raw_mismatch=sum(not r['match'] for r in records),
                   known_native_ub=sum(r['classification'] == 'known_native_ub_empty_map' for r in records),
                   native_crash_safety_divergence=sum(r['classification'] == 'native_crash_safety_divergence' for r in records),
                   actual_fail=sum(r['classification'] == 'actual_fail' for r in records),
                   timeouts=sum(r['classification'] == 'timeout' for r in records),
                   java_exceptions=sum(r['java']['exception'] for r in records),
                   roundtrip_requested=sum(r.get('roundtrip', False) for r in records),
                   roundtrip_executed=sum(r['roundtrip_executed'] for r in records), records=records)
    (results / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf8')
    print(json.dumps({k:v for k,v in summary.items() if k != 'records'}, indent=2))
    for r in records:
        if not r['match']:
            print(r['name'], r['classification'], 'cpp', r['cpp']['rc'], 'java', r['java']['rc'], ','.join(r['diff']))
    return int(summary['actual_fail'] > 0 or summary['timeouts'] > 0
               or summary['java_exceptions'] > 0
               or summary['roundtrip_requested'] != summary['roundtrip_executed'])


if __name__ == '__main__':
    raise SystemExit(main())
