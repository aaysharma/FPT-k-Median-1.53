import json
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent
cert=json.load(open(ROOT/'bfree_final_unified_cert_complete_results.json'))
assert cert['total_cases']==180
assert cert['count']==180
assert cert['bad_count']==0
assert cert['missing']==[]
assert Fraction(cert['min_certificate_value'])>0
assert Fraction(153,100)-Fraction(cert['rho'])==Fraction(1,100000)
manifest=json.load(open(ROOT/'bfree_paper_artifact_manifest.json'))
assert manifest['certificate_status']['failures']==0
print('paper-style artifact summary checks passed')
