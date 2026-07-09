import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
ROOT = Path(__file__).resolve().parent
out = ROOT / 'case179_submit_ready_check.json'
if out.exists():
    out.unlink()
cmd = [sys.executable, 'bfree_final_mixedexp_cert_fast.py', '--start', '179', '--end', '180', '--out', str(out), '--max-depth', '12']
print('$', ' '.join(cmd), flush=True)
subprocess.run(cmd, cwd=ROOT, check=True)
row = json.load(open(out))['rows'][0]
assert row['idx'] == 179
assert row['ok'] is True and row['nfails'] == 0
assert row['nodes'] == 1
assert Fraction(row['min']) == Fraction(2995487733859, 1713034514100000)
assert 'layer-1 GSinf degree-8 Taylor upper' in row.get('M_upper', '')
print('high case 179 verified')
