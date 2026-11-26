**English** | [Русский](REPORT.ru.md)

# Does the Kuznyechik S-box hide an exploitable backdoor? An independent reproduction

Kuznyechik (GOST R 34.12-2015, RFC 7801) uses an 8-bit S-box π, which it shares with
the Streebog hash. The designers say π was chosen at random under quality
criteria. Biryukov–Perrin–Udovenko (2016) and Perrin (2019) showed that π is
instead a highly structured object, a *TKlog*. Stanica–Dutta–Mandal
(arXiv:2507.02181v4) linked that structure to anomalous *extended differentials*
and gave a two-round key-recovery attack.

This report re-checks those claims from scratch and asks one question: **does
the hidden structure give an attacker an advantage that an ordinary S-box
would not?** We did not assume a backdoor exists. Every claim below comes with
code, seeds and raw output in this repository.

## Summary

| Question | Answer |
|---|---|
| Does π have the TKlog structure? | **Yes**, reproduced byte for byte (256/256). |
| Does it cause anomalous differential properties? | **Yes**: extended differential uniformity 64/33/21 vs ≤ 15 for 4000 controls. |
| Does the structure help an attacker beyond equally good S-boxes? | **Only in one narrow sense**: about 4–6 bits in single-characteristic probability at 2–3 rounds (z > 12). Not in the only working attack. |
| Is Kuznyechik weaker than with a random S-box? | **No.** Random S-boxes are easier to attack on every attack metric we measured. |
| Attack on the full 10 rounds? | **None known or found.** The family where the gain appears cannot be used at ≥ 3 rounds at all. |
| Practical key recovery? | Only for **2 rounds**, and that attack works for *any* S-box. |
| Proof of a deliberate backdoor? | **No.** There is strong evidence the published design story is incomplete. That is a different claim. |

## Background: why people suspect a backdoor

1. **Hidden structure.** π can be written as a discrete logarithm in a hidden
   field GF(2⁸)/(X⁸+X⁴+X³+X²+1) plus a few small tables (TKlog). Perrin
   estimates the chance that a random 8-bit permutation is a TKlog at about
   2⁻¹⁶⁰¹.
2. **Conflict with the designers' account.** The designers maintain that π
   came from a pseudo-random search filtered by quality criteria. Perrin reports
   meeting the alleged designer at an ISO meeting in Paris in October 2019. The
   designer kept the randomness claim and said the generating program had been
   lost [who.paris.inria.fr/Leo.Perrin/pi.html](https://who.paris.inria.fr/Leo.Perrin/pi.html). The stated procedure fits the structure that
   was found very poorly. Without the real generator, though, this is strong
   circumstantial evidence and not a proof of deception.
3. **Precedent.** Dual_EC_DRBG showed that unexplained structure in
   agency-supplied constants can hide a trapdoor. That raises the bar for trust;
   it is not evidence of an attack.
4. **Standardisation.** These concerns were raised while Kuznyechik was being
   considered for ISO/IEC 18033-3 (2018–2019). We have not seen the voting
   records, so we do not claim they were the only or the official reason for
   the objections. A 2018 designers' memorandum to WG2 was cited to us but was
   **not verified** here.

## Results

Scripts live in `experiments/verify/` and outputs in `results/verify/`. Run
everything from the repository root. Set `OMP_NUM_THREADS=1` for the parallel
scripts.

### R1. TKlog reconstruction

- **Hypothesis.** Perrin's published TKlog parameters reproduce π exactly.
- **Justification.** π(0) = κ(0). For x = α^l, i = l mod 17 and j = ⌊l/17⌋:
  π(x) = κ(16−j) if i = 0, otherwise κ(16−i) ⊕ (α¹⁷)^{s(j)}.
- **Experiment.** `stage1_2_tklog.py`. Field arithmetic is independent of the
  cipher implementation.
- **Result.** 0 mismatches out of 256. π(0)=0xFC, π(1)=0xEE, and π is a
  permutation. The earlier failure in `experiments/tklog.py` came from the log
  convention. Sage's GF(2⁸) is backed by Givaro, where `_log_repr(1) = 255`, not
  0. With the ordinary log, x=1 gives 0xFC, which collides with π(0). No value
  was patched by hand; the Givaro convention is applied uniformly.
- **Limitations.** None for the reconstruction itself.

### R2. Linear geometry of the TKlog

- **Hypothesis.** The algebraic objects in Perrin's description exist as
  claimed.
- **Experiment.** `stage1_2_tklog.py`. All 256 elements are checked, not
  samples. The per-element table is in `results/verify/tklog_classes.json`.
- **Result.** GF(16) is a subfield and GF(16)* = ⟨α¹⁷⟩. U = span(0x12, 0x26,
  0x24, 0x30) has dimension 4, U ∩ GF(16) = {0} and U ⊕ GF(16) = GF(2⁸). Every
  element falls into its predicted additive layer. The 16 non-trivial
  multiplicative cosets map to 16 distinct additive cosets of GF(16).

### R3. Extended differentials (Stanica–Dutta–Mandal, Tables 1–4)

- **Hypothesis.** N_c(a,b) = #{x : S(c·x ⊕ a) ⊕ S(x) = b}, computed in the
  Kuznyechik field 0x1C3, has the exceptional values reported in the paper.
- **Justification.** Theorem 6 of the paper: the hidden-field transitions
  (120 solutions) carry over with loss at most 2^rank. The rank is the
  F₂-rank of the map x ↦ c·_K x ⊕ c·_T x and equals deg(c). For c = 02 the rank
  is 1, so δ ≥ 60.
- **Experiment.** `stage3_extdiff.py` builds all 255 full tables, one counter
  row per (c,a).
- **Independent checks.**
  - Two field implementations agree on all 65 536 products.
  - A naive count with the cipher's own `gf_mul` agrees on 300 random
    (c,a,b).
  - Every row and column sums to 256 (doubly stochastic).
  - N_c(a,0) = 1 for every c ≠ 1 (Lemma 25).
  - Rows of N_c are pairwise distinct for every c ≠ 1.
  - The identity of Theorem 11 holds for all c.
- **Result.** Every checked number matches the paper:
  - δ_c = 64/33/21 for 02/e1, 04/91 and 03/be; the ordinary DDT gives 8;
  - all of Table 2, including the lower bounds 60/30/8/8/16;
  - the offset splits 64+56, 32+32+32+24 and 21+11;
  - E_c, η_c and I(K;B) in Tables 3 and 4 to 6 decimals;
  - corr(B(c), δ_c) = 0.9788451787 to 10 decimals.

  Two early mismatches were our own error. The paper rounds ⌈h/2^r⌉ up, and we
  had rounded down.

### R4. Exact two-round key recovery (paper §6)

- **Hypothesis.** The two-round attack works on standard Kuznyechik, with
  whitening and the RFC 7801 key schedule.
- **Experiment.** `stage5_two_round.py`.
  - Re-derives σ⊥ = 0.210276, det D*_S ≡ 1 (mod 5), all of Table 5,
    N1 = 14 117, ρ₂ = 0.7074682424211026, N2 = 40 and 240 669 queries.
    Everything matches exactly.
  - Then runs the attack on a vectorised two-round model that is checked
    against the reference implementation.
- **Result.** The full 256-bit master key was recovered for **21/21 keys**: the
  RFC test key plus 20 random keys (seed 2026).
- **Controls.** `control_two_round.py` and `control_large.py --n1` run the same
  attack with π replaced. The best multiplier is taken over all 254 values
  (N1 = 13 884 at c = 4f for π, matching the paper):

  | S-box | N1 (lower = easier to attack) |
  |---|---|
  | π | 13 884 |
  | affine copies A∘π∘B (n=8) | 13 768 / 13 957 / 14 054 (min/median/max) |
  | random permutations (n=8) | 10 999 / 11 191 / 11 347 |

- **Conclusion.** The attack is **generic** for two rounds. Random S-boxes need
  about 20 % less data than π. π sits inside the distribution of its own affine
  copies, so field alignment gives nothing here. The best multiplier, 4f, is
  not one of the structural classes.
- **Limitations.** Only 8+8 controls, because this metric is expensive (about
  27 min on 16 cores). The spread between controls is very narrow.

### R5. Large-scale controls of the S-box metrics

- **Hypothesis.** π is anomalous among S-boxes, and the anomaly comes from
  aligning the hidden field with the linear-layer field.
- **Controls** (`control_large.py --n 1000 --seed 1`):
  - **R**, 1000 uniformly random permutations;
  - **A**, 1000 copies A∘π∘B with independent random invertible affine A and
    B. They have exactly the same DDT and LAT spectra, degree, NL = 100 and
    δ = 8 as π.

  Random S-boxes that meet π's quality criteria (δ = 8) cannot practically be
  sampled; family A is the closest exact substitute. `check_controls.py`
  confirms that A-controls keep the **full TKlog structure** (120 solutions in
  the conjugated field) but not its alignment with the Kuznyechik field.
- **Multiple comparisons.** Each metric is the extreme value over all 254
  multipliers. π and every control are scored the same way, so the comparison
  is max against max.
- **Result.**

  | metric | π | R: min/median/max | A: min/median/max |
  |---|---|---|---|
  | max δ_c | **64** | 9 / 10 / 13 | 9 / 10 / 15 |
  | max E_c | **145 324** | 131 380 / 131 866 / 132 812 | 131 512 / 132 142 / 133 748 |
  | min N2 (stage-2 samples) | 39 | 37 / 38 / 39 | 37 / 38 / 39 |

  π is beyond every one of 2000 controls on the anomaly metrics (p ≤ 1/1001,
  the resolution limit). On attack cost (N2) it is no better than the
  controls, with p = 1.0.

### R6. Characteristic gain over 2–3 rounds (paper §7) against controls

- **Hypothesis.** The extended-differential gain in the 16→1 and 16→1→16
  characteristic families is specific to π.
- **Experiment.** `char_search.py`. It reproduces Tables 6 and 7 exactly:
  c=01 gives −89.1504 / −174.30 and c=02 gives −83.9801 / −169.72. Each S-box is
  scored by its best characteristic over **all** multipliers, c = 1 included.
  The two runs are `--n 500 --seed 3` and `--n 2000 --seed 4`. A separate slow
  pure-Python implementation (`check_controls.py`) agrees to 4 decimals.
- **Result** (seed 4, n = 2000 per family; log₂ probability, higher is better
  for the attacker):

  | | π (c = 02) | R median | A median | π vs A | π vs R |
  |---|---|---|---|---|---|
  | 2 rounds | −83.98 | −85.34 | −89.57 | z = 16.4, 3.9 bits above the best control | 1.6 % of R are better |
  | 3 rounds | −169.72 | −167.30 | −174.89 | z = 12.3, 2.9 bits above the best control | **99.75 % of R are better** |

  Across 4000 controls a non-classical multiplier beat c = 1 only once. The
  seed-3 run gives the same picture.
- **Interpretation.**
  - Against equally good S-boxes, the field alignment gives π a real,
    reproducible gain of about 5–6 bits at the median. The gain is
    deterministic and explained by R3; it is not a statistical fluke.
  - Against random S-boxes, π is still stronger.
- **Why this gain is not an attack.**
  - At 2 rounds a 2⁻⁸⁴ characteristic needs about 2⁸⁴ pairs, while R4 breaks
    2 rounds with 2¹⁸ queries for any S-box.
  - At 3 rounds the best characteristic is 2⁻¹⁶⁹·⁷, below 2⁻¹²⁸. It cannot be
    used even with the full codebook.
  - The paper reports the gain shrinking to about 1.6 bits at 4 rounds.

  This family is therefore closed as a route to key recovery at any number of
  rounds.
- **Limitations.**
  - These are single-characteristic scores in the product model within a
    restricted family.
  - Clustering into differentials, other topologies and a 3-round extension of
    the exact §6 channel remain open. In §6 the second S-box layer already
    shrinks the information by a factor of about 150.

### R7. Nine rounds without pre-whitening (paper §8)

Not reproduced. The authors call it exploratory, because configurations were
selected after preliminary screening. It also omits the initial whitening key,
so it is **not** standard Kuznyechik. It carries no confirmatory weight.

## Issues found in the earlier scripts

- `experiments/sbox_analysis.py` computes nonlinearity from the 8 single output
  bits only. The true NL uses all 255 component functions and is **100**, not 102.
  As a result `results/random_perms.json` compares inflated values.
- The F₂-rank helpers in `experiments/tklog.py` and
  `experiments/extended_differential.py` do not perform correct elimination.
  They happen to return correct ranks on these inputs. `experiments/verify/gf.py`
  has a correct version.
- `tests/` needs `pytest`, which is not installed.

The reference cipher (`src/gost/kuznechik.py`) and the tables in
`data/original/` were not modified.

## Conclusions

1. **TKlog structure: confirmed**, byte for byte, with its full linear
   geometry.
2. **Proven cryptographic weaknesses.** There is an extreme, structurally
   explained anomaly in the extended differential tables (δ_c up to 64). It
   gives about 5 bits of single-characteristic gain at 2–3 rounds compared with
   equally good S-boxes. **No weakness relative to a random S-box was found**:
   on every attack metric we measured, random S-boxes are easier to attack.
3. **Attack on standard 10 rounds: none.** None is published, and none
   followed from this work. The one family where the structure helps is unusable
   from 3 rounds on.
4. **Practical key recovery: only for 2 rounds**, and it is generic. With a
   random S-box it is about 20 % cheaper.
5. **Deliberate backdoor: no evidence.** What *is* established is that π has a
   hidden structure whose chance appearance is negligible, and that the
   designers' account of its origin fits it very poorly. That is a legitimate
   reason for distrust and for continued cryptanalysis. It does not show that
   the structure gives anyone a working trapdoor. The absence of a public attack
   does not show the opposite either: a trapdoor usable only with secret
   knowledge would by design not show up in open analysis.

## Open work

- An upper bound on the information carried by a 3-round extension of the
  exact §6 channel, to see whether any signal survives.
- Differential clustering (beyond single characteristics) for the c = 02/e1
  classes at 3–4 rounds, compared against A-controls.
- A larger N1 control (n ≥ 100 per family; about 4 h on 16 cores).
- An independent audit of this repository by a third party.

## Reproduction

```bash
python3 experiments/verify/stage1_2_tklog.py              # R1, R2
python3 experiments/verify/stage3_extdiff.py              # R3
python3 experiments/verify/stage5_two_round.py 20         # R4 attack, seed 2026
export OMP_NUM_THREADS=1
python3 experiments/verify/control_two_round.py 10        # R4 control
python3 experiments/verify/control_large.py --n 1000 --seed 1        # R5 (~5 min)
python3 experiments/verify/control_large.py --n 8 --seed 98 --n1     # R4 N1 (~27 min)
python3 experiments/verify/char_search.py --n 2000 --seed 4          # R6 (~13 min)
python3 experiments/verify/check_controls.py              # control sanity, slow cross-check
```

## References

- A. Biryukov, L. Perrin, A. Udovenko. *Reverse-Engineering the S-Box of
  Streebog, Kuznyechik and STRIBOBr1.* EUROCRYPT 2016 [ePrint 2016/071](https://eprint.iacr.org/2016/071).
- L. Perrin. *Partitions in the S-Box of Streebog and Kuznyechik.* ToSC 2019(1),
  IACR ePrint 2019/092 [ePrint 2019/092](https://eprint.iacr.org/2019/092).
- L. Perrin. FAQ on the S-box π [who.paris.inria.fr/Leo.Perrin/pi.html](https://who.paris.inria.fr/Leo.Perrin/pi.html).
- P. Stanica, R. Dutta, B. Mandal. *Extended Differential Cryptanalysis of
  Kuznyechik.* arXiv:2507.02181v4 [arXiv:2507.02181](https://arxiv.org/abs/2507.02181).
- [RFC 7801](https://www.rfc-editor.org/rfc/rfc7801), *GOST R 34.12-2015: Block Cipher "Kuznyechik"*.
