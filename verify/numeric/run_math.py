#!/usr/bin/env python3
"""Broad fixed-seed float overload validation, including all boundary bit patterns."""
import collections
import json
import random
import subprocess
import run_numeric as numeric

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
    payload = '\n'.join(line for _, line in all_cases) + '\n'
    (OUT / 'math.cases.txt').write_text(payload)
    cpp = subprocess.run([str(OUT / 'numeric_oracle')], input=payload, text=True, capture_output=True, timeout=120)
    java = subprocess.run(['/usr/bin/java', '-XX:-UsePerfData', '-cp', str(OUT / 'classes'), 'trimal.stat.NumericProbe'],
                          input=payload, text=True, capture_output=True, timeout=120)
    (OUT / 'math.cpp.txt').write_text(cpp.stdout)
    (OUT / 'math.java.txt').write_text(java.stdout)
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
    (OUT / 'math.results.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != 'differences'}, indent=2))
    print('mismatches', len(result['differences']))
    for d in result['differences'][:20]:
        print(json.dumps(d))
    return 1 if result['differences'] or cpp.returncode != 0 or java.returncode != 0 else 0

if __name__ == '__main__':
    raise SystemExit(main())
