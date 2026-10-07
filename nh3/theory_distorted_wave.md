# Why the recoil peak is too small, and the proposed fix: an exact distorted wave for the ejected electron

Theory note (the code is in other files). Notation from O. Zaidi, *Preliminary report on the
molecular extension of the 3CWZ formalism* (CRNA, 2026). Atomic units.

## 1. What the report assumes

The ejected electron moves in the field of the ion, Eq. (11):

$$\Big[-\tfrac12\nabla^2-\frac{Z(r)}{r}\Big]\phi_b=E_b\,\phi_b ,$$

with $Z(r)$ from the spherically averaged Hartree potential, Eqs. (20–21). The 3CWZ model does
not solve this equation. It takes the analytic Coulomb wave, Eq. (9), and puts the local
value $Z(r)$ in the Sommerfeld parameter, Eq. (12): $\eta(r)=-Z(r)/k$. This is exact only
when $Z$ is constant. For NH$_3^+$, $Z_{ion}(r)$ goes from 7 at the N nucleus to 1.8 at
$r=1$ a.u. and 1.01 at $r=4$ a.u. It changes fastest exactly where the orbitals are.

## 2. Test: exact solution of Eq. (11) vs the 3CWZ wave ($E_b=74$ eV, $k=2.33$)

The script is `diagnostics/pw.py`. It integrates the radial equation with the code's own
$Z_{ion}(r)$ (ZSCR).

**Short-range phase shifts (beyond Coulomb $Z=1$):**

| $l$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| $\delta_l$ (rad) | 2.60 | 1.34 | 0.48 | 0.26 | 0.17 | 0.09 | 0.04 |

These are large. They are the shifts that the $Z=1$ wave outside $r=4$ does not contain.

**Partial-wave content.** For a true continuum state, each partial wave $A_l(r)$ has a
**constant** phase $\sigma_l+\delta_l$ (mod $\pi$) at all $r$. Results:
- the exact distorted wave (DW) has a constant phase, as it should;
- the 3CWZ wave has a phase that **drifts with $r$**: for $l=1$, from $-0.96$ to $+2.9$ rad over $r=0.3$–$2.5$. The exact phase is $1.15$ (mod $\pi$);
- the 3CWZ amplitudes are wrong by factors of up to ~7 (s wave at $r=0.3$–$0.6$, p wave at $r=1.5$).

So inside the molecule, where $\langle\phi_b|e^{-i\mathbf p\cdot\mathbf r}|\Phi\rangle$ is
computed, the 3CWZ wave is **not** a solution of Eq. (11).

**Consequence for the recoil peak.** The recoil peak comes from the ejected electron
being back-scattered by the ion. Elastic e$^-$ + NH$_3^+$ cross section from the exact
phase shifts, divided by Rutherford for $Z=1$:

| angle | 30° | 60° | 90° | 120° | 150° | 180° |
|---|---|---|---|---|---|---|
| exact / Rutherford(Z=1) | 2.5 | 4.8 | 4.9 | 29 | 91 | 137 |

The real ion back-scatters about 100 times more than a $Z=1$ Coulomb field. The
variable-charge wave keeps the $Z=1$ phases outside $r\approx4$ and has wrong phases inside, so it
cannot produce this back-scattering. This is consistent with what we found:
- the model recoil is ~10 times too small for every orbital;
- the bad recoil does not come from the numerics (convergence test) or from non-orthogonality
  (orthogonalization made it smaller, see `theory_orthogonalization.md`).

## 3. Proposed fix: distorted wave for $\phi_b$ (keep everything else of 3CWZ)

Solve Eq. (11) exactly with the report's own $Z(r)$, by partial waves:

$$
\phi_b^{-}(\mathbf r)=\sqrt{\tfrac{2}{\pi}}\frac{1}{k r}\sum_{l}(2l+1)\,i^l\,e^{-i(\sigma_l+\delta_l)}\,u_l(r)\,P_l(\hat{\mathbf k}\cdot\hat{\mathbf r}) ,
$$

where $u_l$ is the regular solution normalized to $\sin(kr-l\pi/2+\eta\ln 2kr+\sigma_l+\delta_l)$, with $\eta=-1/k$.

In the code, only the numerical correction of the ejected wave changes. The structure of Eq. (36) stays the same:

$$
\phi_b\Phi\;\to\;\big[\phi_b^{DW}-\phi_b^{C}(Z=1)\big]\Phi\;+\;\phi_b^{C}(Z=1)\Phi .
$$

- The second term is the existing analytic part (TFGA), unchanged.
- The first term replaces QECWF/TECWN. The difference is a finite sum over $l\le l_{max}\approx 10$, since $\delta_l\to0$.
- The phase shifts persist at large $r$, so the difference does **not** vanish beyond $d=4$. Its radial
  grid must cover the orbital, $r\lesssim 12$ a.u. for 3a$_1$, instead of $[0,d]$.
- The projectile and the PCI factor $C_{ab}$ (scattered and incident waves, Fourier transform of $C_{ab}$) are not changed.

**Why this is consistent with the report:**
- It uses exactly the potential of Eqs. (20–21) and the equation of motion, Eq. (11).
- It drops only the approximation $\eta\to\eta(r)$ for the ejected electron.
- It is the ejected-electron part of the distorted-wave approach of Madison and Al-Hagan (report ref. [8]).
- The asymptotic charge is still 1 (Eq. 33), so the Gamow factor and $C_{ab}$ are unchanged.
- The exact wave is (nearly) orthogonal to the bound orbitals, since both see the same
  potential. This also removes the spurious part of $T^0_{dir}$ found in the orthogonalization
  test. The short-range projectile term (`isr=1`) then becomes consistent with the ejected wave.

**Checks planned:**
1. $\delta_l=0$ for all $l$ (pure $Z=1$) must give back the $Z=1$ results.
2. The partial-wave sum of the $Z=1$ Coulomb wave must reproduce the code's `qf1` on the grid.
3. Convergence in $l_{max}$ and in the radial grid.
