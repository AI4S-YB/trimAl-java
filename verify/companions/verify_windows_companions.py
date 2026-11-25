#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Compare Windows Java main entrypoints against Linux native tools via WSL."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'verify' / 'companions'


def windows_java():
    """The Windows JVM reachable from WSL: TRIMAL_WINDOWS_JAVA, or java.exe on PATH."""
    executable = os.environ.get('TRIMAL_WINDOWS_JAVA')
    if executable is None:
        result = subprocess.run(['bash', '-c', 'command -v java.exe'], capture_output=True, text=True)
        executable = result.stdout.strip()
    if not executable:
        raise SystemExit('no Windows java.exe found; set TRIMAL_WINDOWS_JAVA')
    return executable


def run(command):
    try:
        result = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=20)
        return {'rc': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr}
    except subprocess.TimeoutExpired as e:
        return {'rc': 124, 'stdout': e.stdout or b'', 'stderr': e.stderr or b''}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--java-exe', default=windows_java())
    parser.add_argument('--java-jar', type=Path, default=ROOT / 'java' / 'target' / 'trimal.jar')
    parser.add_argument('--out-dir', type=Path,
                        help='where the report is written (default: a new directory under verify/companions/runs/)')
    options = parser.parse_args()
    out_dir = options.out_dir or BASE / 'runs' / time.strftime('%Y%m%d-%H%M%S')
    out_dir.mkdir(parents=True, exist_ok=True)
    java_jar = subprocess.run(['wslpath', '-w', str(options.java_jar.resolve())], check=True,
                              capture_output=True).stdout.decode('utf-8').strip()
    dna = 'verify/companions/inputs/dna.fa'
    aa = 'verify/companions/inputs/aa.fa'
    unaligned = 'verify/companions/inputs/unaligned.fa'
    compare_file = BASE / 'inputs' / 'compareset_windows.txt'
    compare_file.write_bytes(b'verify/companions/inputs/dna.fa\nverify/companions/inputs/dna_alternative.fa\n')
    compare = 'verify/companions/inputs/compareset_windows.txt'
    stats = ['-sgc', '-sgt', '-ssc', '-sst', '-sident', '-scolidentt']
    cases = []
    def add(tool, label, args, output=False):
        cases.append((tool, label, args, output))
    for tool in ['readal', 'statal']:
        add(tool, 'no_args', [])
        add(tool, 'help', ['-h'])
        add(tool, 'version', ['--version'])
        add(tool, 'unknown_rc255', ['--unknown'])
        add(tool, 'missing_value_rc255', ['-in'])
        add(tool, 'missing_file_rc255', ['-in', 'verify/companions/inputs/absent.fa'])
        add(tool, 'empty_file_rc255', ['-in', 'verify/companions/inputs/empty.fa'])
        add(tool, 'no_option_rc255', ['-in', dna])
    for label, flags in [('all_info', ['-format', '-type', '-info']), ('fasta', ['-fasta']),
                         ('reverse', ['-reverse']), ('onlyseqs', ['-onlyseqs']), ('html', ['-html']),
                         ('short_names', ['-phylip_m10'])]:
        add('readal', label, ['-in', dna] + flags)
    for flag in ['-fasta', '-clustal', '-nexus', '-nbrf', '-mega', '-phylip', '-phylip3.2', '-phylip_paml', '-html', '-onlyseqs']:
        add('readal', 'file_' + flag[1:], ['-in', dna, flag], True)
    add('readal', 'unaligned_file_rc255', ['-in', unaligned, '-phylip'], True)
    add('readal', 'mixed_info_file', ['-in', dna, '-fasta', '-format', '-type', '-info'], True)
    add('readal', 'output_directory_rc255', ['-in', dna, '-fasta', '-out', 'verify/companions/inputs'])
    add('statal', 'help_long', ['--help'])
    add('statal', 'version_short', ['-v'])
    add('statal', 'input_alias', ['-i', dna, '-sgc'])
    add('statal', 'dna_all_stats', ['-in', dna] + stats)
    add('statal', 'aa_all_stats', ['-in', aa] + stats)
    add('statal', 'windowed_stats', ['-in', dna, '-w', '1'] + stats)
    add('statal', 'specific_windows', ['-in', aa, '-gw', '1', '-sw', '2'] + stats)
    add('statal', 'unknown_window_rc255', ['-in', dna, '-w', 'text', '-sgc'])
    add('statal', 'single_compare_conflict_rc255', ['-in', dna, '-sfc'])
    add('statal', 'unaligned_rc255', ['-in', unaligned, '-sgc'])
    add('statal', 'all_compare_stats', ['-compareset', compare] + stats + ['-sfc', '-sft'])
    add('statal', 'all_force_stats', ['-forceselect', dna, '-compareset', compare] + stats + ['-sfc', '-sft'])
    add('statal', 'force_window', ['-forceselect', dna, '-compareset', compare, '-cw', '1', '-sfc', '-sft'])
    add('statal', 'custom_matrix', ['-in', dna, '-matrix', 'verify/companions/inputs/identity.matrix', '-ssc', '-sst'])
    results = []
    for tool, label, args, file_output in cases:
        case_dir = BASE / 'windows_outputs' / (tool + '_' + label)
        case_dir.mkdir(parents=True, exist_ok=True)
        destination = case_dir / 'output.aln' if file_output else None
        if destination:
            args = args + ['-out', str(destination.relative_to(ROOT))]
        outputs = {}
        for side in ['native', 'windows_java']:
            if destination and destination.exists():
                destination.unlink()
            if side == 'native':
                command = [str(ROOT / 'verify' / 'bin' / tool)] + args
            else:
                command = [str(options.java_exe), '-XX:-UsePerfData', '-cp', java_jar,
                           'trimal.' + ('ReadAl' if tool == 'readal' else 'StatAl')] + args
            outputs[side] = run(command)
            outputs[side]['file'] = destination.read_bytes() if destination and destination.exists() else None
            for kind in ['stdout', 'stderr', 'file']:
                data = outputs[side][kind]
                if data is not None:
                    (case_dir / (side + '.' + kind)).write_bytes(data)
        diffs = [key for key in ['rc', 'stdout', 'stderr', 'file'] if outputs['native'][key] != outputs['windows_java'][key]]
        row = {'tool': tool, 'name': label, 'args': args, 'passed': not diffs,
               'differences': diffs, 'artifacts': str(case_dir.relative_to(ROOT))}
        for side in outputs:
            row[side] = {'rc': outputs[side]['rc']}
            for kind in ['stdout', 'stderr', 'file']:
                data = outputs[side][kind]
                row[side][kind] = None if data is None else {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
        results.append(row)
        if diffs:
            print('FAIL', tool, label, diffs, flush=True)
    summary = {'total': len(results), 'passed': sum(c['passed'] for c in results),
               'failed': sum(not c['passed'] for c in results)}
    report = {'suite': 'Windows Java CLI entrypoints vs Linux native byte-for-byte',
              'java_exe': str(options.java_exe),
              'java_version': run([str(options.java_exe), '-XX:-UsePerfData', '-version'])['stderr'].decode('utf-8', errors='replace'),
              'java_jar': options.java_jar.resolve().relative_to(ROOT).as_posix(),
              'java_jar_sha256': hashlib.sha256(options.java_jar.read_bytes()).hexdigest(),
              'summary': summary, 'cases': results}
    (out_dir / 'result_windows.json').write_text(json.dumps(report, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    print('REPORT ' + (out_dir / 'result_windows.json').relative_to(ROOT).as_posix(), flush=True)
    return 0 if summary['failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
