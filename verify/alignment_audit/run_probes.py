#!/usr/bin/env python3
"""Direct-source probes for matrix ordering, float stats and locale-neutral parsing."""
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


def main():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.environ['LC_ALL'] = 'C'
    label = sys.argv[1] if len(sys.argv) > 1 else 'round1'
    classes = AUDIT / ('classes_' + label)
    results = AUDIT / ('probes_' + label)
    results.mkdir(parents=True, exist_ok=True)
    cpp = AUDIT / 'audit_probe_cpp'
    objects = ['alignment', 'sequencesMatrix', 'similarityMatrix', 'statisticsConservation', 'statisticsGaps', 'utils']
    commands = [
        ['/usr/bin/g++', '-std=c++11', '-I', str(ROOT / 'upstream/source'), str(AUDIT / 'audit_probe.cpp')]
        + [str(ROOT / ('upstream/source/' + obj + '.o')) for obj in objects] + ['-o', str(cpp)],
        ['/usr/bin/javac', '--release', '11', '-encoding', 'UTF-8', '-cp', str(classes.relative_to(ROOT)), '-d', str(classes), str(AUDIT / 'AuditProbe.java')]
    ]
    for idx, command in enumerate(commands):
        proc = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)
        (results / ('compile%d.stdout' % idx)).write_bytes(proc.stdout)
        (results / ('compile%d.stderr' % idx)).write_bytes(proc.stderr)
        if proc.returncode:
            sys.stdout.buffer.write(proc.stderr)
            raise SystemExit(proc.returncode)
    records = []
    def run(name, cpp_args, java_args=None):
        directory = results / name
        directory.mkdir(exist_ok=True)
        observed = []
        commands = [(['/usr/bin/timeout', '--signal=KILL', '4s', str(cpp)] + cpp_args),
                    ['/usr/bin/java', '-XX:-UsePerfData', '-cp', str(classes.relative_to(ROOT)), 'trimal.AuditProbe'] + (java_args or cpp_args)]
        for impl, command in zip(['cpp', 'java'], commands):
            proc = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=8)
            for stream, data in [('stdout', proc.stdout), ('stderr', proc.stderr)]:
                (directory / (impl + '.' + stream)).write_bytes(data)
            observed.append((proc.returncode, proc.stdout, proc.stderr))
        record = dict(name=name, match=observed[0] == observed[1], cpp_rc=observed[0][0], java_rc=observed[1][0])
        records.append(record)
        return record
    run('matrix_order_names_and_duplicate_rows', ['matrix'])
    locale_file = AUDIT / 'corpus/turkish_locale.nex'
    locale_file.write_bytes(b'#nexus\nbegin data;\ndimensions ntax=2 nchar=4;\nmatrix\nalpha ACGT\nbeta ACGT\n;\nend;\n')
    run('ascii_case_under_turkish_java_locale', ['locale', str(locale_file.relative_to(ROOT))])
    rng = random.Random(1511012)
    files = []
    for idx in range(96):
        n = rng.choice([2, 3, 5, 7, 11, 23, 41])
        length = rng.choice([1, 2, 3, 5, 19, 47, 101, 199])
        file = AUDIT / ('corpus/identity_%03d.fa' % idx)
        file.write_text(''.join('>identity%03d\n%s\n' % (i, ''.join(rng.choice('ACGTNN--acgt') for _ in range(length))) for i in range(n)), encoding='ascii')
        files.append(str(file.relative_to(ROOT)))
    run('identity_96_alignments_precision9', ['identity'] + files)
    contracts = subprocess.run(['/usr/bin/java', '-XX:-UsePerfData', '-cp', str(classes.relative_to(ROOT)), 'trimal.AuditProbe', 'contracts',
                                'verify/alignment_audit/corpus/does_not_exist.fa'], cwd=ROOT,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=8)
    (results / 'contracts.stdout').write_bytes(contracts.stdout)
    (results / 'contracts.stderr').write_bytes(contracts.stderr)
    summary = dict(differential_calls=len(records), identity_alignments=96, match=sum(r['match'] for r in records),
                   contracts_pass=contracts.returncode == 0, records=records)
    (results / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf8')
    print(json.dumps(summary, indent=2))
    return int(summary['match'] != summary['differential_calls'] or not summary['contracts_pass'])


if __name__ == '__main__':
    raise SystemExit(main())
