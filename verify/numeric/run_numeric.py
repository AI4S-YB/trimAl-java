#!/usr/bin/env python3
"""Fixed-seed bit oracle. Run with WSL /usr/bin/python3; no Maven or verify/work."""
import argparse
import collections
import json
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'verify/numeric'
SEED = 15120261005

def fbits(v):
    return struct.unpack('<I', struct.pack('<f', v))[0]

def dbits(v):
    return struct.unpack('<Q', struct.pack('<d', v))[0]

def build():
    cpp = ['utils', 'similarityMatrix', 'statisticsGaps', 'statisticsConservation']
    command = ['g++', '-std=c++11', '-Wall', '-O2', '-I' + str(ROOT / 'upstream/source'),
               str(ROOT / 'verify/oracle/numeric_oracle.cpp')]
    command += [str(ROOT / ('upstream/source/' + f + '.cpp')) for f in cpp]
    subprocess.run(command + ['-o', str(OUT / 'numeric_oracle')], check=True)
    (OUT / 'classes').mkdir(exist_ok=True)
    names = ['util/CppFormat', 'util/CppOut', 'util/Parallel', 'util/Utils', 'util/Values',
             'stat/SimilarityMatrix', 'stat/StatisticsGaps', 'stat/StatisticsConservation']
    subprocess.run(['/usr/bin/javac', '--release', '11', '-encoding', 'UTF-8', '-d', str(OUT / 'classes')]
                   + [str(ROOT / ('java/src/main/java/trimal/' + f + '.java')) for f in names]
                   + [str(OUT / 'NumericProbe.java')], check=True)

def cases():
    rng = random.Random(SEED)
    result = []
    def add(kind, body):
        result.append((kind, f'{kind} {kind}_{sum(k == kind for k, _ in result):04d} {body}'))
    for mode, alphabet in [(1, 'ACGTU'), (3, 'ARNDCQEGHILKMFPSTWYV'), (4, 'ACGTURYKMSWBDHV'),
                           (11, 'ACGTU'), (12, 'ACGTU'), (14, 'ACGTURYKMSWBDHV'), (15, 'ACGTURYKMSWBDHV')]:
        add('M', f'{mode} {alphabet} -')
    matdir = OUT / 'matrices'
    matdir.mkdir(exist_ok=True)
    add('A', '')
    for case in range(180):
        n = rng.randrange(2, 17)
        symbols = ''.join(rng.sample('ABCDEFGHIJKLMNOPQRSTUVWXYZ', n))
        header = ' \t '.join(symbols)
        if case % 7 == 0:
            header = header.lower()
        labelled = case % 2 == 0
        lines = [header]
        order = list(range(n))
        if labelled:
            rng.shuffle(order)
        for row in order:
            vals = [format(rng.uniform(-20, 30) * 10 ** rng.randrange(-3, 4), '.9g') for _ in range(n)]
            if not labelled:
                vals[0] = str(rng.randrange(-100, 101))  # recognizable no-row first token
            lines.append(((symbols[row] + ' ') if labelled else '') + ' '.join(vals))
        ending = '\r\n' if case % 13 == 0 else '\n'
        path = matdir / f'random_{case:04d}.txt'
        path.write_bytes((ending.join(lines) + ending).encode('latin1'))
        add('M', f'0 {symbols} {path}')
    parser_cases = [
        'A C\nA 1 0\nC 0 1\n',
        'a c\n1 0\n0 1\n',
        'A C\r\nA 1 0\r\nC 0 1\r\n',
        'A C\nA 1 0.1234567890123456\nC 0 1\n',
        'A C\nA 1 x\nC 0 1\n',
        'A C\nA 1 2junk\nC 3 4\n',
        'A C\nA 1 1e\nC 3 4\n',
        'A C\nA 1 1e40\nC 3 4\n',
        'A C\nA 1 -1e40\nC 3 4\n',
        'A C\nA 1 0x1p2\nC 3 4\n',
        'A C\nA 1 .25\nC -.25 +2\n',
        'A C\nA 1 NaN\nC 3 4\n',
        'A C\n1 2\na 3\n',
        'A C\n1 2\nA 3\n',
        'A C\nA 1\n',
        'A C\nA 1 2\n',
        'A C\nA 1 2',
        'A A\nA 1 2\nA 3 4\n',
        'A 1\nA 1 2\nC 3 4\n',
    ]
    for case, content in enumerate(parser_cases):
        path = matdir / f'parser_{case:04d}.txt'
        path.write_bytes(content.encode('latin1'))
        add('M', f'0 AC {path}')
    for case in range(700):
        datatype, mode, alphabet = rng.choice([(1, 1, 'ACGTU'), (2, 1, 'ACGU'),
                                               (3, 3, 'ARNDCQEGHILKMFPSTWYV'), (4, 4, 'ACGTURYKMSWBDHV'),
                                               (5, 14, 'ACGTURYKMSWBDHV')])
        n, cols = rng.randrange(1, 25), rng.randrange(1, 121)
        hw = rng.randrange(cols // 4 + 1)
        gap_prob = rng.uniform(0, 1)
        indet = 'X' if datatype == 3 else 'N'
        rows = []
        for _ in range(n):
            row = ''
            for _ in range(cols):
                p = rng.random()
                row += '-' if p < gap_prob else (indet if rng.random() < .12 else rng.choice(alphabet))
            rows.append(row)
        minimum = rng.choice([0, 25, 33.3333333, 50, 60, 99.9, 100])
        sim = rng.choice([0, .00001, .01, .1, .5, .9, 1])
        gap = rng.choice([0, .01, .1, .2, .5, .9, 1])
        hw2 = rng.randrange(cols // 4 + 1) if case % 5 == 0 else -1
        add('S', f'{datatype} {n} {cols} {hw} {mode} {fbits(minimum)} {fbits(sim)} {fbits(gap)} {hw2} ' + ' '.join(rows))
    for case in range(1400):
        type_ = rng.choice(['f', 'd'])
        width = 32 if type_ == 'f' else 64
        raw = rng.getrandbits(width)
        p = rng.choice([0, 1, 2, 4, 6, 10, 17])
        fixed = rng.randrange(2)
        add('F', f'{type_} {raw} {p} {fixed}')
    for type_, raws in [('f', [0, 0x80000000, 0x7fc00000, 0xffc00000, 0xffa01234, 0x7fa01234, fbits(1.4142135), fbits(math.exp(-.25))]),
                        ('d', [0, 1, 0x8000000000000000, 0x7ff8000000000000, 0xfff8000000000000, 0x7ff0000000000001, 0xfff0000000000001])]:
        for raw in raws:
            for p in [0, 1, 6, 10, -1]:
                for fixed in [0, 1]:
                    add('F', f'{type_} {raw} {p} {fixed}')
    parse = ['', ' ', '\t\r\n\v\f-23tail', '\r123', '\v123', '\xa0123', '\x1c12', '000001234567890',
             '2147483647', '2147483648', '2147483649', '4294967295', '4294967296', '-2147483649',
             '9223372036854775807', '9223372036854775808', '-9223372036854775809', '9' * 200,
             '+1.5e3tail', '-.5e+3x', '.e4', '+', '-e', '2e', '2e+', '1.2.3', '1e400', '-1e400',
             'inf', 'infinity!', '-INFINITY', 'NaN', '-nan(payload)', '+nan(abc)', '0x1p2', '-0x1.8p+2',
             '0x.8p-4', '0x1', '0x1p', '-0', '1\x00e5', '12\x00x', '12e2\r', '.000001', '++++', '--2']
    parse += [rng.choice(['', ' ', '\r', '\v', '\t']) + rng.choice(['', '+', '-'])
              + str(rng.getrandbits(rng.randrange(1, 110))) + rng.choice(['', 'tail', '.5', 'e3']) for _ in range(200)]
    for s in parse:
        # Empty C strings can be represented by NUL in the batch protocol.
        add('I', (s if s else '\0').encode('latin1').hex())
    for v in [-3.5, -2.5, -1.5, -.5, 0, .5, 1.5, 2.5, 2147483647., 2147483647.5, 2147483648.,
              -2147483648., -2147483649., math.inf, -math.inf, math.nan] + [rng.uniform(-10000, 10000) for _ in range(200)]:
        add('R', str(dbits(v)))
    for case in range(300):
        type_ = rng.choice(['f', 'i', 'p'])
        n = rng.randrange(0, 90)
        if type_ == 'f':
            vals = [fbits(rng.choice([math.nan, -0., 0., rng.uniform(-10, 10)])) for _ in range(n)]
        else:
            vals = [rng.choice([1, 1, 3, rng.randrange(-30, 30), 16777216, 16777217, 16777218]) for _ in range(n)]
        add('Q', f'{type_} {n} ' + ' '.join(map(str, vals)))
    for case in range(5000):
        v = rng.uniform(-40, 40) if case % 2 == 0 else 10 ** rng.uniform(-30, 30)
        add('T', str(fbits(v)))
    return result

def parse_lines(output):
    return {line.split('|')[0]: dict(item.split('=', 1) for item in line.split('|')[1:])
            for line in output.splitlines()}

def compare(all_cases, native, java):
    totals, matches, diffs = collections.Counter(), collections.Counter(), []
    for kind, line in all_cases:
        id_ = line.split()[1]
        totals[kind] += 1
        left, right = native.get(id_), java.get(id_)
        # Double results remain in the raw output as diagnostics. T is testing
        # the float overload: double libm vs Math is not its equality gate.
        if kind == 'T':
            left = {k: v for k, v in (left or {}).items() if k not in ('expDouble', 'log10Double')}
            right = {k: v for k, v in (right or {}).items() if k not in ('expDouble', 'log10Double')}
        if left == right:
            matches[kind] += 1
        else:
            fields = []
            if left is not None and right is not None:
                for field in set(left) | set(right):
                    if left.get(field) != right.get(field):
                        fields.append(field)
            diffs.append({'id': id_, 'kind': kind, 'input': line, 'fields': sorted(fields), 'cpp': left, 'java': right})
    return {'seed': SEED, 'counts': dict(totals), 'matches': dict(matches), 'total': len(all_cases),
            'exact_matches': sum(matches.values()), 'differences': diffs}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-build', action='store_true')
    parser.add_argument('--label', default='results')
    args = parser.parse_args()
    if not args.no_build:
        build()
    all_cases = cases()
    payload = '\n'.join(line for _, line in all_cases) + '\n'
    (OUT / 'cases.txt').write_text(payload)
    cpp = subprocess.run([str(OUT / 'numeric_oracle')], input=payload, text=True, capture_output=True, timeout=120)
    java = subprocess.run(['/usr/bin/java', '-XX:-UsePerfData', '-cp', str(OUT / 'classes'), 'trimal.stat.NumericProbe'],
                          input=payload, text=True, capture_output=True, timeout=120)
    (OUT / f'{args.label}.cpp.txt').write_text(cpp.stdout)
    (OUT / f'{args.label}.java.txt').write_text(java.stdout)
    (OUT / f'{args.label}.cpp.stderr.txt').write_text(cpp.stderr)
    (OUT / f'{args.label}.java.stderr.txt').write_text(java.stderr)
    result = compare(all_cases, parse_lines(cpp.stdout), parse_lines(java.stdout))
    result['exit_codes'] = {'cpp': cpp.returncode, 'java': java.returncode}
    native_rows, java_rows = parse_lines(cpp.stdout), parse_lines(java.stdout)
    diagnostic = collections.Counter()
    for kind, line in all_cases:
        if kind != 'T': continue
        id_ = line.split()[1]
        left, right = native_rows.get(id_, {}), java_rows.get(id_, {})
        for field in ('expDouble', 'log10Double'):
            if left.get(field) != right.get(field): diagnostic[field] += 1
    result['diagnostic_double_differences_not_gates'] = dict(diagnostic)
    result['production_mismatches'] = len(result['differences'])
    result['skipped_cases'] = 0
    (OUT / f'{args.label}.json').write_text(json.dumps(result, indent=2))
    summary = {k: v for k, v in result.items() if k != 'differences'}
    summary['difference_fields'] = dict(collections.Counter(f for d in result['differences'] for f in d['fields']))
    print(json.dumps(summary, indent=2))
    return 1 if result['differences'] or cpp.returncode != 0 or java.returncode != 0 else 0

if __name__ == '__main__':
    sys.exit(main())
