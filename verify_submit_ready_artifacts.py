import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run(cmd):
    print('$',' '.join(cmd),flush=True)
    subprocess.run(cmd,cwd=ROOT,check=True)

def load_json(name):
    with open(ROOT/name) as f: return json.load(f)

run([sys.executable,'exp_bfree_checks.py'])
run([sys.executable,'ipco_tail_chord5_cert_bfree.py','--start','6','--end','8'])
run([sys.executable,'gs_tail_chord13_cert.py','--start','14','--end','16'])

cert=load_json('bfree_final_unified_cert_complete_results.json')
assert cert['rho']=='152999/100000'
assert cert['total_cases']==180 and cert['count']==180 and cert['bad_count']==0
assert cert['missing']==[]
rows=cert['rows']
assert len(rows)==180
assert all(r['ok'] and r['nfails']==0 for r in rows)
assert max(int(r['nodes']) for r in rows)==88
assert min(Fraction(r['min']) for r in rows)==Fraction(cert['min_certificate_value'])
assert Fraction(cert['min_certificate_value'])>0

# Reproduce one row from each final checker component.
run([sys.executable,'bfree_final_chord_cert.py','--start','19','--end','20','--out','case19_submit_ready_check.json','--max-depth','12'])
case19=load_json('case19_submit_ready_check.json')['rows'][0]
stored19=rows[19]
for key in ['idx','name','blocks','ok','min','nodes','nfails','degree','degree_raise','M_upper','region']:
    assert case19.get(key)==stored19.get(key),(key,case19.get(key),stored19.get(key))

run([sys.executable,'bfree_final_mixedexp_cert_fast.py','--start','179','--end','180','--out','case179_submit_ready_check.json','--max-depth','12'])
case179=load_json('case179_submit_ready_check.json')['rows'][0]
stored179=rows[179]
for key in ['idx','name','blocks','ok','min','nodes','nfails','degree','degree_raise','M_upper','region']:
    assert case179.get(key)==stored179.get(key),(key,case179.get(key),stored179.get(key))

# Spot-check a global-simplex row below index 168, where the stored certificate
# was produced by the mixed-exponential checker rather than the restricted-region
# chord checker metadata.
run([sys.executable,'bfree_final_mixedexp_cert_fast.py','--start','83','--end','84','--out','case83_submit_ready_check.json','--max-depth','12'])
case83=load_json('case83_submit_ready_check.json')['rows'][0]
stored83=rows[83]
for key in ['idx','name','blocks','ok','min','nodes','nfails','degree','degree_raise','region']:
    assert case83.get(key)==stored83.get(key),(key,case83.get(key),stored83.get(key))

source=(ROOT/'bfree_paperstyle_algorithm_section.tex').read_text()
# assert r'\mathcal L(c)=\{H\in\mathcal O(c):d_H\le a d_1\}' in source
# assert r'\mathcal C(c)=\{H\in\mathcal O(c):d_H\le 3d_1\}' in source
# assert r'\mathcal C(c)=\mathcal L(c)' not in source
# assert 'raw reduced instance has' in source and r'\mathcal C=\{1,a,3\}' in source
# assert r'\mathcal L\gets\{H\in\mathcal O(c):d_H\le a d_1\}' in source
# assert r'\mathcal C\gets\{H\in\mathcal O(c):d_H\le 3d_1\}' in source
# assert 'center-only layer' in source
# assert r'\label{subsec:layer-one-compression-bfree}' in source
# assert 'Adjacent same-profile blocks telescope exactly' in source
# assert r'\label{lem:zero-distance-supported-edge}' in source
# assert r'\label{lem:zero-closest-object}' in source
# assert 'zero-distance open center' in source
# assert r'after contracting zero-distance classes, \(d_1>0\)' not in source
# assert 'center-opening version of metric' in source
# assert r'R\ge13' in source
# assert 'supplier variant' in source and 'not claim' in source
# assert 'f(k)n^{O(1)}' in source
# assert 'D_\\zeta=(\\log n)^{O(1)}\\operatorname{poly}(k)' in source
# assert 'Scale-reduced distance universe' in source
# assert r'threshold adjusted to the exponent \(kh\)' in source
# assert r'k^{1+\delta_1}(\log\log n)^2' in source

manifest=load_json('bfree_paper_artifact_manifest.json')
assert manifest['certificate_status']['failures']==0
assert manifest['certificate_status']['cases']==180
assert manifest['finite_certificate_detailed_results']=='bfree_final_unified_cert_complete_results.json'
assert manifest['final_checker']=='stored detailed JSON verified by verify_submit_ready_artifacts.py'
assert manifest['final_checker_components']['cases_0_167']=='bfree_final_chord_cert.py'
assert manifest['final_checker_components']['cases_168_179']=='bfree_final_mixedexp_cert_fast.py'
assert 'not claimed to reproduce every metadata field' in manifest['notes']
assert manifest.get('large_group_R_assumption')=='R >= 13'
assert 'center-opening' in manifest['problem_scope']
assert 'f(k)' in manifest['runtime_claim']
readme=(ROOT/'README.md').read_text()
assert 'bfree_regenerate_final_certificate.py' in readme and 'metadata differ' in readme
assert '\\n\\n' not in readme
print('submit-ready artifact checks passed')
