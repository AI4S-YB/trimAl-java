#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Broad fixed-seed float overload validation, including all boundary bit patterns."""
import collections
import json
import random
import subprocess
import run_numeric as numeric
import toolchain

OUT = numeric.OUT
SEED = numeric.SEED + 1

def main():
    numeric.build()
    rng = random.Random(SEED)
    raws = [rng.getrandbits(32) for _ in range(120000)]
    raws += [numeric.fbits(rng.uniform(-104, 89)) for _ in range(60000)]
    raws += [numeric.fbits(-rng.uniform(0, 35)) for _ in range(60000)]
    raws += [numeric.fbits(10 ** rng.uniform(-44, 38)) for _ in range(40000)]
    for center in [0, 1, 0x007fffff, 0x00800000, 0x3f800000, numeric.fbits(88.0),
                   0x7f800000, 0x7fc00000, 0xff800000, 0xffc00000,
                   numeric.fbits(float.fromhex('0x1.62e42ep6')),
                   numeric.fbits(-float.fromhex('0x1.9fe368p6'))]:
        raws += [(center + i) & 0xffffffff for i in range(-512, 513)]
    all_cases = [('T', f'T T_{i:06d} {raw}') for i, raw in enumerate(raws)]

    # double pow() at the port's two call sites. Site 1 (the strict method's similarity
    # cut) receives the windowed conservation values of two columns, i.e. floats in
    # [0, 1]; site 2 (the standard deviation of the column identities) receives a
    # column identity k / sequences and the float average of such identities.
    units = [0.0, 1.0, 1e-45, 1.17549435e-38, 1.1920929e-07, 0.5, 0.99999994, 0.1, 0.2, 0.25, 0.33333334]
    pairs = []
    for i in range(100000):
        a = rng.choice(units) if i % 4 == 0 else rng.uniform(0, 1)
        b = rng.choice(units) if i % 8 == 1 else rng.uniform(0, 1)
        pairs.append(('W1', f'W W1_{i:06d} 1 {numeric.fbits(a)} {numeric.fbits(b)}'))
    for i in range(100000):
        n, k = rng.randrange(1, 200), None
        identity = rng.randrange(0, n + 1) / n
        avg = (rng.randrange(0, n + 1) / n) if i % 2 else rng.uniform(0, 1)
        pairs.append(('W2', f'W W2_{i:06d} 2 {numeric.fbits(identity)} {numeric.fbits(avg)}'))
    all_cases += pairs
    payload = '\n'.join(line for _, line in all_cases) + '\n'
    run = numeric.run_dir('math')
    (run / 'cases.txt').write_text(payload)
    cpp = subprocess.run([str(OUT / 'numeric_oracle')], input=payload, text=True, capture_output=True, timeout=120)
    java = subprocess.run([toolchain.java_tool('java'), '-XX:-UsePerfData', '-cp', str(OUT / 'classes'), 'trimal.stat.NumericProbe'],
                          input=payload, text=True, capture_output=True, timeout=120)
    (run / 'cpp.txt').write_text(cpp.stdout)
    (run / 'java.txt').write_text(java.stdout)
    result = numeric.compare(all_cases, numeric.parse_lines(cpp.stdout), numeric.parse_lines(java.stdout))
    result['seed'] = SEED
    result['exit_codes'] = {'cpp': cpp.returncode, 'java': java.returncode}
    cpp_rows, java_rows = numeric.parse_lines(cpp.stdout), numeric.parse_lines(java.stdout)
    diagnostic = collections.Counter()
    for kind, line in all_cases:
        id_ = line.split()[1]
        c, j = cpp_rows.get(id_, {}), java_rows.get(id_, {})
        for field in ['exp', 'log10']:
            if c.get(field) != c.get(field + 'Double'):
                diagnostic['cpp_float_vs_double_' + field] += 1
            if c.get(field + 'Double') != j.get(field + 'Double'):
                diagnostic['cpp_vs_java_double_' + field] += 1
    result['diagnostic_differences_not_gates'] = dict(diagnostic)
    result['production_mismatches'] = len(result['differences'])
    result['skipped_cases'] = 0
    (run / 'results.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != 'differences'}, indent=2))
    print('REPORT ' + str(run / 'results.json'))
    print('mismatches', len(result['differences']))
    for d in result['differences'][:20]:
        print(json.dumps(d))
    return 1 if result['differences'] or cpp.returncode != 0 or java.returncode != 0 else 0

if __name__ == '__main__':
    raise SystemExit(main())
