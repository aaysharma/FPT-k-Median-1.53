#!/usr/bin/env python3
"""Regenerate a valid current-format 180-case Bernstein certificate.

The packaged detailed JSON is an authoritative certificate artifact and may
contain historical metadata such as a global-simplex region label. This wrapper
regenerates the same mathematical inequalities in the current checker format:
  * cases 0..167: bfree_final_chord_cert.py
  * cases 168..179: bfree_final_mixedexp_cert_fast.py

The latter is required because these high-index rows contain a layer-one
GS2^infty@1 block and use a degree-eight Taylor upper bound for that layer-one
exponential, matching the paper text. The output is intended to be a valid
certificate, not a byte-for-byte reproduction of every stored metadata field.
"""
import argparse, json, os, subprocess, sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run(cmd):
    print('$',' '.join(cmd),flush=True)
    subprocess.run(cmd,cwd=ROOT,check=True)

def load(p):
    with open(p) as f: return json.load(f)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='regenerated_180_cases.json')
    ap.add_argument('--max-depth',default='12')
    ap.add_argument('--keep-parts',action='store_true')
    args=ap.parse_args()
    p1=ROOT/'_regen_0_168_chord.json'
    p2=ROOT/'_regen_168_180_mixedexp.json'
    run([sys.executable,'bfree_final_chord_cert.py','--start','0','--end','168','--out',str(p1),'--max-depth',str(args.max_depth)])
    run([sys.executable,'bfree_final_mixedexp_cert_fast.py','--start','168','--end','180','--out',str(p2),'--max-depth',str(args.max_depth)])
    rows=load(p1)['rows']+load(p2)['rows']
    rows.sort(key=lambda r:r['idx'])
    ids={r['idx'] for r in rows}
    bad=[r for r in rows if (not r.get('ok')) or r.get('nfails',0)!=0]
    out={
        'rho':'152999/100000',
        'description':'Final B-free reduced 180-case certificate; indices 0-167 use layer-a chord checker, indices 168-179 use layer-a chord plus layer-1 GSinf degree-8 Taylor upper.',
        'total_cases':180,
        'count':len(rows),
        'bad_count':len(bad),
        'missing':[i for i in range(180) if i not in ids],
        'max_nodes':max(int(r.get('nodes',0)) for r in rows),
        'min_certificate_value':str(min(Fraction(r['min']) for r in rows)),
        'rows':rows,
        'bad':bad,
    }
    with open(ROOT/args.out,'w') as f: json.dump(out,f,indent=2)
    print(json.dumps({k:out[k] for k in ['count','bad_count','missing','max_nodes','min_certificate_value']},indent=2))
    if not args.keep_parts:
        for p in (p1,p2):
            try: os.remove(p)
            except FileNotFoundError: pass
if __name__=='__main__': main()
