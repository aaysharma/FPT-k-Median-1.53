import subprocess
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def run(cmd):
    print('$', ' '.join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)

run([sys.executable, 'ipco_tail_chord5_cert_bfree.py', '--start', '6', '--end', '8'])
run([sys.executable, 'gs_tail_chord13_cert.py', '--start', '14', '--end', '16'])
print('tail spotchecks passed')
