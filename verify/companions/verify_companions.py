#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Dedicated readAl/statAl byte-for-byte differential suite, run inside WSL."""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'verify'))
BASE = ROOT / 'verify' / 'companions'
INPUTS = BASE / 'inputs'
OUTPUTS = BASE / 'outputs'
RUNS = BASE / 'runs'
CLASSES = BASE / 'classes'
JAR = BASE / 'companions.jar'
NATIVE = ROOT / 'verify' / 'bin'
STATS = ['-sgc', '-sgt', '-ssc', '-sst', '-sident', '-scolidentt']
FORMATS = ['-clustal', '-fasta', '-fasta_m10', '-nbrf', '-nexus', '-mega',
           '-phylip3.2', '-phylip3.2_m10', '-phylip', '-phylip_m10',
           '-phylip_paml', '-phylip_paml_m10', '-html', '-onlyseqs']

import toolchain  # noqa: E402


def repo_path(path):
    return os.path.relpath(Path(path).resolve(), ROOT)


def write_fasta(name, sequences, names=None):
    path = INPUTS / (name + '.fa')
    if names is None:
        names = ['alpha_long_name', 'beta_long_name', 'gamma_long_name', 'delta_long_name'][:len(sequences)]
    path.write_bytes(''.join('>' + n + ' description ' + str(i + 1) + '\n' + s + '\n'
                             for i, (n, s) in enumerate(zip(names, sequences))).encode('ascii'))
    return repo_path(path)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def difference(a, b):
    if a == b:
        return None
    offset = next((i for i in range(min(len(a), len(b))) if a[i] != b[i]), min(len(a), len(b)))
    return {'offset': offset, 'native_length': len(a), 'java_length': len(b),
            'native_context': repr(a[max(0, offset - 40):offset + 100]),
            'java_context': repr(b[max(0, offset - 40):offset + 100])}


def run_process(command):
    try:
        result = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                timeout=20, env=dict(os.environ, LC_ALL='C', TZ='UTC'))
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as e:
        return 124, e.stdout or b'', e.stderr or b''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-compile', action='store_true')
    parser.add_argument('--java-jar', type=Path, help='Use an existing integrated jar and skip compilation')
    parser.add_argument('--filter', default='')
    parser.add_argument('--out-dir', type=Path,
                        help='where the report is written (default: a new directory under verify/companions/runs/)')
    options = parser.parse_args()
    java_jar = options.java_jar.resolve() if options.java_jar else JAR
    out_dir = options.out_dir or RUNS / time.strftime('%Y%m%d-%H%M%S')
    out_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    INPUTS.mkdir(parents=True, exist_ok=True)
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    if not options.no_compile and options.java_jar is None:
        if CLASSES.exists():
            shutil.rmtree(CLASSES)
        CLASSES.mkdir(parents=True)
        sources = sorted((ROOT / 'java' / 'src' / 'main' / 'java').rglob('*.java'))
        result = subprocess.run([toolchain.java_tool('javac'), '--release', '11', '-encoding', 'UTF-8', '-d', str(CLASSES)]
                                + [str(p) for p in sources], cwd=ROOT, capture_output=True)
        (BASE / 'compile.stdout').write_bytes(result.stdout)
        (BASE / 'compile.stderr').write_bytes(result.stderr)
        if result.returncode:
            sys.stderr.buffer.write(result.stderr)
            return result.returncode
        subprocess.run([toolchain.tool('jar'), '--create', '--file', str(JAR), '-C', str(CLASSES), '.'], check=True, cwd=ROOT)

    dna = write_fasta('dna', ['ACTGACTG--ACGT--', 'ACT-A-TGGTACG-T-', 'ACTGAC-G-TACGT--'])
    aa = write_fasta('aa', ['ARNDCQEGHILKMFPSTWYV--', 'AR-DCEEGH-LKMFP-TWYVAA', 'ARNDCQEGHIL-MF-STWY---'])
    rna = write_fasta('rna', ['ACGUACGU--ACGU--', 'ACU-A-UGGUACG-U-', 'ACUGAC-G-UACGU--'])
    deg = write_fasta('degenerate_dna', ['ACGTRYSWKMBDHVN-', 'ACGT-N-ACGTACG--', 'ACGTACGTACGTACGT'])
    rdeg = write_fasta('degenerate_rna', ['ACGURYSWKMBDHVN-', 'ACGU-N-ACGUACG--', 'ACGUACGUACGUACGU'])
    unaligned = write_fasta('unaligned', ['ACTG--ACT', 'ACT', 'ACTGACTGACT'])
    single = write_fasta('single', ['ACTG--ACTG'])
    indet = write_fasta('indeterminate', ['NNNNNNNNNNNN', 'NNN---NNNNNN', 'N-NNNNNNN-NN'])
    gaps = write_fasta('gaps', ['---ACTG---ACT---', '---ACTGA--ACT---', '---ACTGG--AC----'])
    short = write_fasta('short', ['AT', 'AC', 'AG'])
    empty = INPUTS / 'empty.fa'; empty.write_bytes(b'')
    invalid = INPUTS / 'invalid.txt'; invalid.write_bytes(b'not an alignment\n')
    absent = repo_path(INPUTS / 'absent.fa')
    matrix = INPUTS / 'identity.matrix'
    matrix.write_bytes(b'A C G T U\nA 1 0 0 0 0\nC 0 1 0 0 0\nG 0 0 1 0 0\nT 0 0 0 1 1\nU 0 0 0 1 1\n')
    invalid_matrix = INPUTS / 'invalid.matrix'; invalid_matrix.write_bytes(b'A A\n1 0\n0 1\n')
    missing_symbols = INPUTS / 'missing_symbols.matrix'; missing_symbols.write_bytes(b'A C\n1 0\n0 1\n')
    comp2 = write_fasta('dna_alternative', ['ACTG-ACTGACGT---', 'ACTATG-GTAC-GT--', 'ACTGACGTACGT----'])
    comp3 = write_fasta('dna_alternative_order', ['ACTGACGT--ACGT--', 'ACTGACTGACGT----', 'ACTATGGTACGT----'],
                        ['gamma_long_name', 'alpha_long_name', 'beta_long_name'])
    comp_paths = INPUTS / 'compareset.txt'
    comp_paths.write_bytes(('\n'.join([dna, comp2, comp3]) + '\n').encode('ascii'))
    set_empty = INPUTS / 'compareset_empty.txt'; set_empty.write_bytes(b'')
    set_bad = INPUTS / 'compareset_bad.txt'; set_bad.write_bytes((dna + '\n' + absent + '\n').encode('ascii'))
    set_unaligned = INPUTS / 'compareset_unaligned.txt'; set_unaligned.write_bytes((dna + '\n' + unaligned + '\n').encode('ascii'))
    set_type = INPUTS / 'compareset_type.txt'; set_type.write_bytes((dna + '\n' + aa + '\n').encode('ascii'))
    set_single = INPUTS / 'compareset_single.txt'; set_single.write_bytes((dna + '\n').encode('ascii'))

    cases = []
    def add(tool, name, args, file_output=False):
        args = [repo_path(value) if isinstance(value, Path)
                or (isinstance(value, str) and value.startswith(str(ROOT))) else value
                for value in args]
        cases.append({'tool': tool, 'name': name, 'args': args, 'file_output': file_output})
    for tool in ['readal', 'statal']:
        for name, args in [('no_arguments', []), ('help', ['-h']), ('version', ['--version']),
                           ('unknown', ['--unknown']), ('missing_in_value', ['-in']),
                           ('help_with_extra', ['-h', '-in', dna]), ('missing_file', ['-in', absent]),
                           ('empty_file', ['-in', str(empty)]), ('invalid_file', ['-in', str(invalid)]),
                           ('in_only', ['-in', dna]), ('duplicate_in', ['-in', dna, '-in', dna])]:
            add(tool, name, args)
    add('readal', 'help_long_invalid', ['--help'])
    add('readal', 'version_short_invalid', ['-v'])
    for name, data in [('dna', dna), ('aa', aa), ('rna', rna), ('deg', deg), ('rdeg', rdeg),
                       ('unaligned', unaligned), ('single', single), ('gaps', gaps), ('short', short), ('indet', indet)]:
        for flag in ['-format', '-type', '-info']:
            add('readal', name + flag, ['-in', data, flag])
        add('readal', name + '_info_combination', ['-in', data, '-format', '-type', '-info'])
    for flag in FORMATS:
        add('readal', 'stdout_' + flag[1:], ['-in', dna, flag])
        add('readal', 'file_' + flag[1:], ['-in', dna, flag], True)
        add('readal', 'reverse_' + flag[1:], ['-in', dna, flag, '-reverse'])
    for name, data in [('aa', aa), ('rna', rna), ('deg', deg), ('unaligned', unaligned)]:
        for flag in ['-fasta', '-nbrf', '-onlyseqs']:
            add('readal', name + '_' + flag[1:], ['-in', data, flag])
    add('readal', 'reverse_default', ['-in', dna, '-reverse'])
    add('readal', 'reverse_info_output', ['-in', dna, '-reverse', '-format', '-type', '-info'], True)
    add('readal', 'all_info_output', ['-in', dna, '-fasta', '-format', '-type', '-info'], True)
    for name, args in [('missing_input', ['-fasta']), ('missing_output_value', ['-in', dna, '-out']),
                       ('duplicate_format', ['-in', dna, '-fasta', '-nexus']),
                       ('duplicate_info', ['-in', dna, '-info', '-info']),
                       ('duplicate_reverse', ['-in', dna, '-reverse', '-reverse']),
                       ('info_stdout_conflict', ['-in', dna, '-info', '-fasta']),
                       ('reverse_info_conflict', ['-in', dna, '-info', '-reverse']),
                       ('info_output_conflict', ['-in', dna, '-info', '-out', str(OUTPUTS / 'unused')]),
                       ('output_no_option', ['-in', dna, '-out', str(OUTPUTS / 'unused')]),
                       ('output_directory', ['-in', dna, '-fasta', '-out', str(INPUTS)]),
                       ('output_missing_parent', ['-in', dna, '-fasta', '-out', str(BASE / 'absent_dir' / 'out.fa')])]:
        add('readal', name, args)
    for flag in ['-clustal', '-nexus', '-phylip', '-phylip3.2', '-phylip_paml', '-mega']:
        add('readal', 'unaligned_file_' + flag[1:], ['-in', unaligned, flag], True)

    # Native-generated format fixtures also exercise Java's existing format loaders.
    for flag in ['-clustal', '-nbrf', '-nexus', '-mega', '-phylip', '-phylip3.2', '-phylip_paml']:
        target = INPUTS / ('native_' + flag[1:] + '.aln')
        rc, out, err = run_process([repo_path(NATIVE / 'readal'), '-in', dna, flag, '-out', repo_path(target)])
        if rc != 0:
            raise RuntimeError('Unable to prepare ' + str(target) + ': ' + repr(err))
        add('readal', 'input_' + flag[1:], ['-in', str(target), '-format', '-type', '-info'])
        add('readal', 'roundtrip_' + flag[1:], ['-in', str(target), flag])
        add('statal', 'input_' + flag[1:], ['-in', str(target)] + STATS)

    add('statal', 'help_long', ['--help'])
    add('statal', 'version_short', ['-v'])
    add('statal', 'input_alias', ['-i', dna, '-sgc'])
    for name, data in [('dna', dna), ('aa', aa), ('rna', rna), ('deg', deg), ('rdeg', rdeg),
                       ('gaps', gaps), ('short', short)]:
        for flag in STATS:
            add('statal', name + '_' + flag[1:], ['-in', data, flag])
        add('statal', name + '_all_stats', ['-in', data] + STATS)
        add('statal', name + '_out_of_order_stats', ['-in', data] + list(reversed(STATS)))
    for name, args in [('duplicates_stats', ['-in', dna, '-sgc', '-sgc']),
                       ('unaligned_stats', ['-in', unaligned, '-sgc']),
                       ('stats_missing_input', ['-sgc']), ('sfc_missing_set', ['-sfc']),
                       ('sft_single_input', ['-in', dna, '-sft']),
                       ('single_and_set', ['-in', dna, '-compareset', str(comp_paths), '-sgc']),
                       ('single_and_reference', ['-in', dna, '-forceselect', dna, '-sgc']),
                       ('two_incompatibilities', ['-in', dna, '-compareset', str(comp_paths), '-sfc']),
                       ('force_missing_set', ['-forceselect', dna, '-sgc']),
                       ('general_and_specific_window', ['-in', dna, '-w', '1', '-gw', '1', '-sgc']),
                       ('matrix_missing', ['-in', dna, '-matrix', absent, '-ssc']),
                       ('matrix_invalid', ['-in', dna, '-matrix', str(invalid_matrix), '-ssc']),
                       ('matrix_missing_symbols', ['-in', dna, '-matrix', str(missing_symbols), '-ssc']),
                       ('alternative_invalid', ['-in', dna, '--alternative_matrix', 'unknown', '-ssc']),
                       ('set_missing', ['-compareset', absent, '-sfc']),
                       ('set_empty', ['-compareset', str(set_empty), '-sfc']),
                       ('set_bad_input', ['-compareset', str(set_bad), '-sfc']),
                       ('set_unaligned', ['-compareset', str(set_unaligned), '-sfc']),
                       ('set_datatype_mismatch', ['-compareset', str(set_type), '-sfc'])]:
        add('statal', name, args)
    for flag in ['-w', '-gw', '-sw', '-cw']:
        for value in ['0', '-1', 'text', '1.9', '1e2', '100', '']:
            add('statal', 'window_' + flag[1:] + '_' + (value or 'empty'), ['-in', dna, flag, value] + STATS)
        add('statal', 'window_' + flag[1:] + '_missing', ['-in', dna, flag])
        add('statal', 'window_' + flag[1:] + '_duplicate', ['-in', dna, flag, '1', flag, '2', '-sgc'])
        add('statal', 'window_' + flag[1:] + '_valid', ['-in', dna, flag, '1'] + STATS)
    add('statal', 'specific_windows', ['-in', aa, '-gw', '1', '-sw', '2'] + STATS)
    add('statal', 'custom_matrix', ['-in', dna, '-matrix', str(matrix), '-ssc', '-sst'])
    add('statal', 'custom_matrix_no_similarity', ['-in', dna, '-matrix', str(matrix), '-sgc'])
    add('statal', 'custom_and_alternative', ['-in', dna, '-matrix', str(matrix), '--alternative_matrix',
                                           'degenerated_nt_identity', '-ssc', '-sst'])
    for name, data in [('dna', dna), ('deg', deg), ('rdeg', rdeg), ('aa', aa)]:
        add('statal', 'alternative_' + name, ['-in', data, '--alternative_matrix',
                                             'degenerated_nt_identity', '-ssc', '-sst'])
    for label, compare_set in [('three', comp_paths), ('one', set_single)]:
        for flag in ['-sfc', '-sft']:
            add('statal', 'compare_' + label + '_' + flag[1:], ['-compareset', str(compare_set), flag])
        add('statal', 'compare_' + label + '_all', ['-compareset', str(compare_set)] + STATS + ['-sfc', '-sft'])
    for flag in ['-w', '-cw']:
        add('statal', 'compare_window_' + flag[1:], ['-compareset', str(comp_paths), flag, '1', '-sfc', '-sft'])
    add('statal', 'compare_no_stats', ['-compareset', str(comp_paths)])
    add('statal', 'force_all_stats', ['-forceselect', dna, '-compareset', str(comp_paths)] + STATS + ['-sfc', '-sft'])
    add('statal', 'force_window', ['-forceselect', dna, '-compareset', str(comp_paths), '-cw', '1', '-sfc', '-sft'])
    add('statal', 'force_reversed_args', ['-compareset', str(comp_paths), '-forceselect', dna, '-sfc', '-sft'])

    # Thirty independent small datasets, with deterministic seeds and mixed alphabets.
    rng = random.Random(151)
    for index in range(30):
        alphabet = ['ACGT', 'ACGU', 'ARNDCQEGHILKMFPSTWYV'][index % 3]
        count = 3 + index % 4
        length = 8 + index % 17
        seqs = [''.join(rng.choice(alphabet + '--') for _ in range(length)) for _ in range(count)]
        data = write_fasta('random_%02d' % index, seqs, ['seq_%d' % n for n in range(count)])
        add('statal', 'random_%02d_all' % index, ['-in', data] + STATS)
        add('statal', 'random_%02d_windows' % index, ['-in', data, '-w', '1'] + STATS)
        add('readal', 'random_%02d_info' % index, ['-in', data, '-format', '-type', '-info'])
        add('readal', 'random_%02d_reverse' % index, ['-in', data, '-reverse', '-fasta'])

    # Byte-transparent headers, CRLF input, missing final LF, mixed case, and
    # indeterminate/gap punctuation exercise the existing shared alignment layer.
    byte_fixtures = {
        'crlf': (ROOT / dna).read_bytes().replace(b'\n', b'\r\n'),
        'no_final_lf': (ROOT / dna).read_bytes().rstrip(b'\n'),
        'utf8_headers': ('>s\u00e9quence_\u00e1\u00e9\u00ed_long \u4e2d\u6587\u63cf\u8ff0\nACTGACTG--ACGT--\n'
                         '>\u03b2eta_long_\u6d4b\u8bd5 description\nACT-A-TGGTACG-T-\n'
                         '>third_\u65e5\u672c\u8a9e_name desc\nACTGAC-G-TACGT--\n').encode('utf-8'),
        'latin1_headers': b'>s\xe9q_long_name description\nACTGACTG--ACGT--\n>other_\xff_name desc\nACT-A-TGGTACG-T-\n>third_name desc\nACTGAC-G-TACGT--\n',
        'lowercase': b'>alpha\nactgactg--acgt--\n>beta\naCt-A-tGgtacG-t-\n>gamma\nACTGAC-g-TACGT--\n',
        'punctuation_nt': b'>alpha\nACTG?N..--ACGT--\n>beta\nACT-?N..GTACG-T-\n>gamma\nACTGAN..-TACGT--\n',
        'punctuation_aa': b'>alpha\nARNDCX??GHILMFPSTWYV--\n>beta\nAR-DCX??G-ILMFP-TWYVAA\n>gamma\nARNDCX??GHILMF-STWY---\n',
        'all_gap_column': b'>alpha\nACTG--ACTG--ACTG\n>beta\nACTG--AC-G--ACTG\n>gamma\nACTG--AC-G--AC-G\n',
        'all_gap_row': b'>alpha\nACTG--ACTG--ACTG\n>beta\nACTG--AC-G--ACTG\n>gamma\n----------------\n',
    }
    for label, content in byte_fixtures.items():
        fixture = INPUTS / ('bytes_' + label + '.fa')
        fixture.write_bytes(content)
        add('readal', 'bytes_' + label + '_info', ['-in', str(fixture), '-format', '-type', '-info'])
        add('readal', 'bytes_' + label + '_fasta', ['-in', str(fixture), '-fasta'])
        add('readal', 'bytes_' + label + '_m10_file', ['-in', str(fixture), '-phylip_m10'], True)
        add('readal', 'bytes_' + label + '_html', ['-in', str(fixture), '-html'], True)
        add('statal', 'bytes_' + label + '_all_stats', ['-in', str(fixture)] + STATS)
        add('statal', 'bytes_' + label + '_window', ['-in', str(fixture), '-w', '1'] + STATS)
    for label, content in [('crlf', comp_paths.read_bytes().replace(b'\n', b'\r\n')),
                           ('no_final_lf', comp_paths.read_bytes().rstrip(b'\n'))]:
        fixture = INPUTS / ('compareset_' + label + '.txt')
        fixture.write_bytes(content)
        add('statal', 'bytes_compareset_' + label, ['-compareset', str(fixture), '-sfc', '-sft'])
    add('statal', 'bytes_relative_path', ['-in', dna] + STATS)
    add('readal', 'bytes_relative_path', ['-in', dna, '-format', '-type', '-info'])

    if options.filter:
        cases = [c for c in cases if options.filter in c['tool'] + ':' + c['name']]
    results = []
    for index, case in enumerate(cases):
        case_dir = OUTPUTS / (case['tool'] + '_' + case['name'].replace('.', '_'))
        case_dir.mkdir(parents=True, exist_ok=True)
        args = case['args'][:]
        destination = None
        if case['file_output']:
            destination = case_dir / 'output.aln'
            args += ['-out', repo_path(destination)]
        native_command = [repo_path(NATIVE / case['tool'])] + args
        java_command = [toolchain.java_tool('java'), '-XX:-UsePerfData', '-cp', repo_path(java_jar), 'trimal.' + ('ReadAl' if case['tool'] == 'readal' else 'StatAl')] + args
        side_results = {}
        for side, command in [('native', native_command), ('java', java_command)]:
            if destination and destination.exists():
                destination.unlink()
            rc, stdout, stderr = run_process(command)
            file_bytes = destination.read_bytes() if destination and destination.exists() else None
            (case_dir / (side + '.stdout')).write_bytes(stdout)
            (case_dir / (side + '.stderr')).write_bytes(stderr)
            if file_bytes is not None:
                (case_dir / (side + '.aln')).write_bytes(file_bytes)
            side_results[side] = {'rc': rc, 'stdout': stdout, 'stderr': stderr, 'file': file_bytes}
        native, java = side_results['native'], side_results['java']
        diffs = {key: difference(native[key], java[key]) for key in ['stdout', 'stderr']}
        diffs['rc'] = None if native['rc'] == java['rc'] else {'native': native['rc'], 'java': java['rc']}
        if native['file'] is None or java['file'] is None:
            diffs['file'] = None if native['file'] == java['file'] else {'native_exists': native['file'] is not None, 'java_exists': java['file'] is not None}
        else:
            diffs['file'] = difference(native['file'], java['file'])
        passed = all(value is None for value in diffs.values())
        # Do not call a divergence PASS. Native statAl deletes uninitialized
        # pointers after the first CRLF list path fails (statAl.cpp:365/468).
        # Classify only the exact observed SIGABRT and identical emitted error;
        # any changed output, timeout, Java success, or other failure remains FAIL.
        native_abort = b'double free or corruption (out)\n'
        known_native_ub = (
            case['tool'] == 'statal' and case['name'] == 'bytes_compareset_crlf'
            and native['rc'] == -6 and java['rc'] == 255
            and native['stdout'] == java['stdout'] and native['file'] == java['file']
            and native['stderr'] == java['stderr'] + native_abort
            and java['stderr'].startswith(b'\nERROR: Alignment not loaded: "')
            and b'dna.fa\r". Check input file content\n\n' in java['stderr'])
        status = 'pass' if passed else 'known_native_ub' if known_native_ub else 'fail'
        row = {'tool': case['tool'], 'name': case['name'], 'args': args, 'passed': passed,
               'status': status,
               'differences': {key: value for key, value in diffs.items() if value is not None},
               'artifacts': repo_path(case_dir), 'native': {}, 'java': {}}
        for side in ['native', 'java']:
            row[side] = {'rc': side_results[side]['rc']}
            for key in ['stdout', 'stderr', 'file']:
                data = side_results[side][key]
                row[side][key] = None if data is None else {'bytes': len(data), 'sha256': digest(data)}
        results.append(row)
        if not passed:
            print(status.upper(), case['tool'], case['name'], json.dumps(row['differences'], ensure_ascii=True), flush=True)
        elif (index + 1) % 50 == 0:
            print('Processed', index + 1, '/', len(cases), flush=True)
    summary = {'total': len(results), 'passed': sum(c['passed'] for c in results),
               'failed': sum(c['status'] == 'fail' for c in results),
               'nonmatching': sum(not c['passed'] for c in results),
               'known_divergence': sum(c['status'] == 'known_native_ub' for c in results),
               'by_tool': {tool: {'total': sum(c['tool'] == tool for c in results),
                                  'passed': sum(c['tool'] == tool and c['passed'] for c in results),
                                  'failed': sum(c['tool'] == tool and c['status'] == 'fail' for c in results),
                                  'known_divergence': sum(c['tool'] == tool and c['status'] == 'known_native_ub' for c in results)}
                           for tool in ['readal', 'statal']},
               'elapsed_seconds': round(time.time() - start, 2)}
    report = {'suite': 'readAl/statAl v1.5.1 byte-for-byte differential', 'summary': summary,
              'java_version': subprocess.run([toolchain.java_tool('java'), '-XX:-UsePerfData', '-version'], capture_output=True).stderr.decode('utf-8'),
              'compile': None if options.no_compile or options.java_jar is not None else
                         'javac --release 11 -encoding UTF-8 -d verify/companions/classes <all current Java sources>',
              'native_version': {tool: run_process([repo_path(NATIVE / tool), '--version'])[1].decode('latin1') for tool in ['readal', 'statal']},
              'java_jar': repo_path(java_jar), 'java_jar_sha256': digest(java_jar.read_bytes()),
              'seed': 151, 'independent_random_datasets': 30, 'cases': results}
    report_name = 'result.json' if not options.filter else 'result_' + options.filter.replace(':', '_') + '.json'
    (out_dir / report_name).write_text(json.dumps(report, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    print('REPORT ' + repo_path(out_dir / report_name), flush=True)
    return 0 if summary['failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
