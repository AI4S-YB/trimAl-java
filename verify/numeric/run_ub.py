#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Observe UB separately: original output is never silently counted as a crash/skip."""
import collections
import hashlib
import json
import os
import subprocess
import run_numeric as numeric
import toolchain

OUT = numeric.OUT
ROOT = numeric.ROOT

def invoke(binary, payload, env=None):
    try:
        r = subprocess.run(binary, input=payload, text=True, capture_output=True, timeout=3, env=env)
        return {'rc': r.returncode, 'stdout': r.stdout, 'stderr': r.stderr,
                'stdout_sha256': hashlib.sha256(r.stdout.encode()).hexdigest()}
    except subprocess.TimeoutExpired as e:
        return {'timeout': True, 'stdout': str(e.stdout), 'stderr': str(e.stderr)}

def main():
    numeric.build()
    cpp = [str(OUT / 'numeric_oracle')]
    java = [toolchain.java_tool('java'), '-XX:-UsePerfData', '-cp', str(OUT / 'classes'), 'trimal.stat.NumericProbe']
    result = {'seed': numeric.SEED, 'aa_alternative': {}, 'all_threshold': {}, 'empty_gaps': {}}
    payload = 'A alternative_aa\n'
    natives = [invoke(cpp, payload) for _ in range(12)]
    j = invoke(java, payload)
    result['aa_alternative'] = {'native_runs': natives, 'java': j,
                               'all_native_outputs_identical': len({r.get('stdout_sha256') for r in natives}) == 1,
                               'exact_matches': sum(r.get('stdout') == j.get('stdout') and r.get('rc') == j.get('rc') for r in natives)}
    # Every value >= threshold makes calcCutPoint read vectAux[-1].
    cases = []
    for cols in [1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 100, 1000, 33000, 50000]:
        for value, threshold in [(1, 0), (1, 1), (.5, .5), (.25, .2)]:
            cases.append(f'U U_{cols}_{value}_{threshold} {cols} {numeric.fbits(value)} {numeric.fbits(threshold)}')
    payload = '\n'.join(cases) + '\n'
    natives = [invoke(cpp, payload) for _ in range(6)]
    j = invoke(java, payload)
    result['all_threshold'] = {'cases_per_run': len(cases), 'native_runs': natives, 'java': j,
                              'all_native_outputs_identical': len({r.get('stdout_sha256') for r in natives}) == 1,
                              'exact_matches': sum(r.get('stdout') == j.get('stdout') and r.get('rc') == j.get('rc') for r in natives)}
    gap_records = []
    for n in [0, 1, 2, 5, 20]:
        payload = f'G empty_{n} {n}\n'
        gap_records.append({'species': n, 'native_runs': [invoke(cpp, payload) for _ in range(3)], 'java': invoke(java, payload)})
    result['empty_gaps'] = gap_records
    sanitizer = OUT / 'numeric_oracle_asan'
    command = ['g++', '-std=c++11', '-Wall', '-O2', '-g', '-fsanitize=address,undefined',
               '-fno-sanitize-recover=all', '-I' + str(ROOT / 'upstream/source'),
               str(ROOT / 'verify/oracle/numeric_oracle.cpp')]
    command += [str(ROOT / ('upstream/source/' + name + '.cpp')) for name in
                ['utils', 'similarityMatrix', 'statisticsGaps', 'statisticsConservation']]
    compile_result = subprocess.run(command + ['-o', str(sanitizer)], text=True, capture_output=True)
    result['sanitizer_build'] = {'rc': compile_result.returncode, 'stderr': compile_result.stderr}
    if compile_result.returncode == 0:
        env = dict(os.environ, ASAN_OPTIONS='detect_leaks=0')
        result['sanitizers'] = {
            'aa_alternative': invoke([str(sanitizer)], 'A alt_asan\n', env),
            'all_threshold': invoke([str(sanitizer)], f'U cut_asan 8 {numeric.fbits(1)} 0\n', env),
            'empty_gaps': invoke([str(sanitizer)], 'G gaps_asan 1\n', env)}
    result['counts'] = {
        'aa_original_invocations': 12, 'aa_exact_java_matches': result['aa_alternative']['exact_matches'],
        'cut_original_invocations': 6, 'cut_cases': len(cases) * 6,
        'cut_run_exact_java_matches': result['all_threshold']['exact_matches'],
        'empty_gap_original_invocations': 15,
        'empty_gap_output_runs': sum(bool(r.get('stdout')) for g in gap_records for r in g['native_runs']),
        'empty_gap_crash_runs': sum(r.get('rc', 0) < 0 for g in gap_records for r in g['native_runs']),
        'empty_gap_timeouts': sum(r.get('timeout', False) for g in gap_records for r in g['native_runs'])}
    result['known_ub_output_differences'] = [
        {'kind': 'empty_gaps', 'species': g['species'], 'native_stdout': r.get('stdout'),
         'java_stdout': g['java'].get('stdout'), 'native_rc': r.get('rc'), 'java_rc': g['java'].get('rc')}
        for g in gap_records for r in g['native_runs']
        if r.get('stdout') != g['java'].get('stdout') or r.get('rc') != g['java'].get('rc')]
    result['counts']['known_ub_output_differences'] = len(result['known_ub_output_differences'])
    result['counts']['sanitizer_confirmed_ub_categories'] = sum(
        r.get('rc', 0) != 0 for r in result.get('sanitizers', {}).values())
    run = numeric.run_dir('ub')
    (run / 'results.json').write_text(json.dumps(result, indent=2))
    (run / 'cases.txt').write_text('\n'.join(cases) + '\n')
    print(json.dumps(result['counts'], indent=2))
    print('Sanitizers', {k: v.get('rc') for k, v in result.get('sanitizers', {}).items()})
    print('REPORT ' + str(run / 'results.json'))
    # Empty-gap outputs are intentionally recorded as UB, not matched/ignored
    # defined cases. Gate the measured AA layout, fallback and sanitizer evidence.
    failed = result['aa_alternative']['exact_matches'] != 12 or result['all_threshold']['exact_matches'] != 6
    failed |= result['counts']['sanitizer_confirmed_ub_categories'] != 3
    failed |= any(g['java'].get('rc') != 0 for g in gap_records)
    return 1 if failed else 0

if __name__ == '__main__':
    raise SystemExit(main())
