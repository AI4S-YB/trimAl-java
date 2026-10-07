#!/usr/bin/env python3
"""Run owned JUnit tests in verify/numeric/classes using cached Maven jars."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import run_numeric as numeric


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--m2-repository', type=Path, default=os.environ.get('MAVEN_REPOSITORY'))
    options = parser.parse_args()
    numeric.build()
    artifacts = [('org/junit/jupiter', 'junit-jupiter-api', '5.10.2'),
                 ('org/junit/jupiter', 'junit-jupiter-engine', '5.10.2'),
                 ('org/junit/platform', 'junit-platform-commons', '1.10.2'),
                 ('org/junit/platform', 'junit-platform-engine', '1.10.2'),
                 ('org/junit/platform', 'junit-platform-launcher', '1.10.2'),
                 ('org/apiguardian', 'apiguardian-api', '1.1.2'),
                 ('org/opentest4j', 'opentest4j', '1.3.0')]
    candidates = ([options.m2_repository] if options.m2_repository else
                  [Path.home() / '.m2/repository'] + list(Path('/mnt/c/Users').glob('*/.m2/repository')))
    cache = next((candidate for candidate in candidates if all(
        (candidate / group / artifact / version / f'{artifact}-{version}.jar').is_file()
        for group, artifact, version in artifacts)), None)
    if cache is None:
        raise SystemExit('JUnit dependencies not cached; run mvn -f java/pom.xml test first '
                         'or pass --m2-repository PATH')
    jars = [cache / group / artifact / version / f'{artifact}-{version}.jar' for group, artifact, version in artifacts]
    classpath = ':'.join([str(numeric.OUT / 'classes')] + list(map(str, jars)))
    tests = ['util/CppFormatTest', 'util/UtilsQuicksortTest', 'util/UtilsTest',
             'stat/SimilarityMatrixTest', 'stat/StatisticsGapsTest', 'stat/StatisticsConservationTest']
    subprocess.run(['/usr/bin/javac', '--release', '11', '-encoding', 'UTF-8', '-cp', classpath,
                    '-d', str(numeric.OUT / 'classes')]
                   + [str(numeric.ROOT / ('java/src/test/java/trimal/' + t + '.java')) for t in tests]
                   + [str(numeric.OUT / 'NumericJUnit.java')], check=True)
    run = subprocess.run(['/usr/bin/java', '-XX:-UsePerfData', '-cp', classpath, 'NumericJUnit'], capture_output=True, text=True)
    (numeric.OUT / 'junit.results.txt').write_text(run.stdout + run.stderr)
    summary = {'exit_code': run.returncode, 'tests': tests}
    for name, count in re.findall(r'\[\s*(\d+) tests (found|started|successful|failed)\s*\]', run.stdout):
        summary[count] = int(name)
    (numeric.OUT / 'junit.results.json').write_text(json.dumps(summary, indent=2))
    print(run.stdout + run.stderr)
    return run.returncode

if __name__ == '__main__':
    raise SystemExit(main())
