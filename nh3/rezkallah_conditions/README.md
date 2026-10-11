# Tests of the NH₃ BBK3CW-Z code under the conditions of Rezkallah et al., NIM B 269, 2750 (2011)

The experimental data and the paper's 1CW/PWBA curves were read directly from the paper's
vector figures (`data/`, files `EMS_Bawagan_*`, `EMS_Rezkallah1CW_*`, `DDCS_Opal_*`,
`DDCS_Rezkallah1CW_*`). All model runs use `nh3_opt_claude.for` with only the kinematics
changed (`scripts/mk.py`, `scripts/mkd.py`). The 1e results include both components
(1ex + 1ey, factor 2).

## 1. EMS (Bawagan et al. 1988): symmetric non-coplanar, E₀ = 1200 eV + IP, 600/600 eV, θ = 45°

Run with `igeom=2`: scattered electron at θs = 45°, φs = 180°; ejected electron at θe = 45°;
φe = Φ from 0 to 30° in 1° steps. Figure: `EMS_NH3_1200eV_symmetric_noncoplanar_45deg.png`.

- **Peak positions** are correct for all three orbitals: q ≈ 0.58 (1e), 0.1 (2a₁) and 0.57 (3a₁) a.u.
- **Absolute value** (default grids) is about 0.57 times the paper's PWBA for all three orbitals.
  Part of this is numerical: with finer grids (nleg = 16, md = 4) the 3a₁ peak goes from
  0.74 to 1.00 (×10⁻³ a.u.), against 1.32 for the PWBA. The default grids are **not converged**
  when the ejected electron is fast (600 eV, k = 6.6 a.u.). The changes go up to −40%/+70% along the profile.
- **Shape** (model × 1.75):
  - 2a₁ and 3a₁ follow the data up to q ≈ 2.4.
  - 1e is about 1.5 times below the data for q > 1.3. The high-momentum tail of the
    Moccia 1e orbital is too small, as in the paper's own PWBA curve.
  - At q < 0.2 the data for the p-like 1e and 3a₁ are 2 to 7 times above every calculation
    (a known effect of distortion and finite resolution in EMS).

## 2. DDCS (Opal et al. 1972): Ei = 500 eV, θe = 120°, versus ejection energy

The TDCS is integrated over the scattered electron direction:
- in θs, 17 points equally spaced in ln K, Simpson rule;
- in φ, 0–180° in 45° steps, Simpson rule.

The three valence orbitals are summed. Figure: `DDCS_NH3_Ei500eV_thetae120deg.png`;
data: `data/DDCS_BBK3CWZ_NH3_Ei500_the120.dat`.

| Ee (eV) | model / Opal | share from θs < ~20–35° |
|---|---|---|
| 5 | 1.79 | 0.82–0.90 |
| 20 | 1.58 | 0.60–0.89 |
| 60 | 1.35 | 0.34–0.64 |
| 150 | 7.1 | 0.26–0.31 |

- **Ee = 5–60 eV**: the model is 1.3–1.8 times above the data. The paper's 1CW curve is
  closer to the data.
- **Ee = 150 eV**: the excess comes from large scattering angles (large momentum transfer K),
  where the TDCS with the default grids levels off instead of falling. It is **numerical**:
  - finer grids (nleg = 16, nd = 8, md = 4) lower the TDCS at θs = 68° by a factor of 12–30;
  - at θs = 120° even the finer grid is not stable;
  - orthogonalizing the ejected wave (`iorth=1`) does not remove it (6.4 instead of 7.1).
- In the region that dominates at low Ee (small θs), finer grids change the TDCS by ±25%.

## Conclusions

1. The Moccia orbitals and the orientation average are fine near the momentum-profile peak.
   The 1e high-momentum tail is too small (factor ~1.5 at q > 1.3).
2. The numerical part of the code (TSCWN/TECWN on r < d) is converged only when the outgoing
   electrons are slow and K is small. For fast electrons (600 eV) or large K (> 3 a.u.), the
   default grids `nleg=10, nd=4, md=1` give errors from tens of percent up to a factor of 30.
   The El Mir kinematics (θs = 6°, K ≈ 0.85, Ee = 74 eV) are in the safe region (±20%), so the
   earlier conclusions about the recoil peak stand.
3. Integrated quantities (DDCS, SDCS, TCS) need grids adapted to K, or the analytic
   (1CW-like) limit at large K, before they can be compared absolutely.
