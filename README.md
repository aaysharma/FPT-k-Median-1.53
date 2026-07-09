# Submit-ready unified B-free proof package

Problem scope: centers may be opened at client points; the supplier variant where client points are forbidden as centers is not claimed.

This package contains a paper-style LaTeX/PDF writeup and the computational artifacts for the unified B-free soft-LP local factor certificate for the center-opening version of general-metric k-median (client points may be opened as centers).  The final assignment rule uses leaders up to a d_1 and direct centers up to 3 d_1; the center-only layer is analytically reduced to the layer-a envelope before the finite certificate.

Assignment rule used in the paper:

```text
L(c) = {G : d_G <= (3/(2+eps)) d_1}
C(c) = {G : d_G <= 3 d_1}
```

The analysis eliminates the escape cap, reduces raw distances to layers `1,a,3`, pushes the center-only `3` layer down to `a`, and dominates all layer-`a` material by one `GS2^infty@a` envelope.


Main files:

- `bfree_paperstyle_full.pdf`: compiled paper-style PDF.
- `bfree_paperstyle_full.tex`: wrapper source.
- `bfree_paperstyle_algorithm_section.tex`: main paper section with algorithms, structural lemmas, assignment analysis, B-free transformation, and final certificate theorem.
- `bfree_final_unified_cert_complete_results.json`: detailed row-by-row result file for the final 180-case Bernstein certificate.
- `bfree_final_unified_cert_complete_results.json`: authoritative detailed row-by-row certificate artifact.
- `bfree_regenerate_final_certificate.py`: optional regeneration wrapper for a valid 180-case certificate in the current checker format. It runs `bfree_final_chord_cert.py` for cases 0--167 and `bfree_final_mixedexp_cert_fast.py` for cases 168--179. It may produce semantically equivalent rows whose metadata differ from the packaged historical JSON.
- `bfree_final_chord_cert.py`: checker component for cases 0--167, where the layer-a chord bound suffices.
- `bfree_final_mixedexp_cert_fast.py`: checker component for cases 168--179, where a layer-one `GSinf@1` block is present and a degree-eight Taylor bound is used for that layer-one exponential.
- `ipco_tail_chord5_cert_bfree.py`: IPCO multiplicity/tail chord verifier.  Supports `--start` and `--end` ranges.
- `gs_tail_chord13_cert.py`: GS2 multiplicity/tail chord verifier.  Supports `--start` and `--end` ranges.
- `exp_bfree_checks.py`: exponential and residual tail margin checks.
- `verify_submit_ready_artifacts.py`: fast package consistency verifier.
- `verify_high_case179.py`: optional regeneration of high closure case `idx=179` using the final mixed-exponential checker.

Fast verification:

```bash
python3 verify_submit_ready_artifacts.py
```

This audits the detailed JSON certificate, reruns the exponential checks, spot-checks the tail scripts, checks the corrected assignment rule/source scope, and audits representative certificate rows from both checker components. To regenerate the memory-heavy high closure case `idx=179`, run `python3 verify_high_case179.py`.

Optional regeneration of a valid final 180-case certificate in the current checker format:

```bash
python3 bfree_regenerate_final_certificate.py --max-depth 12 --out regenerated_180_cases.json
```

The packaged `bfree_final_unified_cert_complete_results.json` is the authoritative detailed certificate artifact. The optional regeneration command is intended to reproduce the same mathematical inequalities, but row metadata such as region labels or node counts may differ from the stored historical certificate.

Full optional tail checks:

```bash
python3 ipco_tail_chord5_cert_bfree.py --start 6 --end 60
python3 gs_tail_chord13_cert.py --start 14 --end 60
```

The paper text uses the detailed `bfree_final_unified_cert_complete_results.json` file.  This archive no longer relies only on a summary artifact.


Zero-distance convention:
The paper does not rely on zero-distance contraction to imply d_1>0.  Instead, Lemmas `Zero-distance supported edge` and `Zero closest ordinary object` handle the degenerate case: if the closest compressed ordinary object has averaged distance zero, a zero-distance center is open in every rounding outcome, so the client is omitted from the normalized factor program.
