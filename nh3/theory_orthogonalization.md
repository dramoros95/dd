# Orthogonalization of the ejected-electron wave in the molecular 3CWZ (BBK3CW-Z) model

Theory note for the option `iorth=1` of `nh3_opt_claude.for` (the code is in that file;
this file contains only the theory). Notation follows O. Zaidi, *Preliminary report on
the molecular extension of the 3CWZ formalism* (CRNA, 2026). Atomic units.

## 1. Why orthogonalization is needed

The direct amplitude is (report, Eqs. 23–26)

$$
T_{dir}=T^1_{dir}-T^0_{dir},\qquad
T^0_{dir}=\Big\langle \phi_a^-\,\phi_b^-\,C_{ab}^-\Big|\frac{Z(r_0)}{r_0}\Big|\phi_i^+\,\Phi\Big\rangle ,
$$

where $Z(r_0)=1$ in the frozen-core form of Eq. (7) and $Z(r_0)=Z_{ion}(r_0)$ when the
short-range part of the projectile–core potential is kept (option `isr=1`).

If the ejected-electron wave $\phi_b$ were an eigenfunction of the same one-electron
Hamiltonian as the bound orbital $\Phi$, we would have $\langle\phi_b|\Phi\rangle=0$.
In first order (no $C_{ab}$), $T^0_{dir}$ would then vanish identically, because the
operator $Z(r_0)/r_0$ does not act on the target electron:

$$
T^0_{dir}\big|_{C_{ab}=1}=\langle\phi_a|Z/r_0|\phi_i\rangle\,\langle\phi_b|\Phi\rangle .
$$

In the 3CWZ model, $\phi_b$ is a Coulomb wave with the variable charge $Z_b(r)$, while
$\Phi$ is a Moccia SCF orbital. They are **not** orthogonal, so $T^0_{dir}$ contains a
spurious contribution proportional to $\langle\phi_b|\Phi\rangle$. With `isr=1`, this
spurious part is multiplied by $Z_{ion}(r)-1$, which is large (up to 6) exactly where
the s-like density of the orbital is concentrated. The NH$_3$ runs show this:
- the recoil peak becomes too large for 3a$_1$ and 2a$_1$;
- the binary peak of 2a$_1$ is distorted.

The standard cure is to orthogonalize the ejected wave to the bound state (Schmidt):

$$
\tilde\phi_b=\phi_b-\Phi\,\frac{\langle\Phi|\phi_b\rangle}{\langle\Phi|\Phi\rangle}.
$$

## 2. Where the ejected wave enters the code

After the Fourier decomposition of $C_{ab}$ (report, Eqs. 30–32), the ejected electron
enters only through

$$
F(\mathbf p)=(2\pi)^{-3/2}\int \phi_b^{-*}(\mathbf r)\,e^{-i\mathbf p\cdot\mathbf r}\,\chi(\mathbf r)\,d^3r ,
$$

with $\chi$ the ionized orbital (or one of its channels, §3). In the code,
$F(\mathbf p)=$ `tfv(icm)` (analytic $Z=1$ part, TFGA) $+$ `qsm(icm,k)`
(numerical TECWN part, $r<d$).

Replacing $\phi_b$ by $\tilde\phi_b$ gives

$$
\tilde F(\mathbf p)=F(\mathbf p)-\frac{S\;G(\mathbf p)}{n},\qquad
S=\langle\phi_b|\chi\rangle=(2\pi)^{3/2}F(\mathbf 0),
$$

$$
G(\mathbf p)=(2\pi)^{-3/2}\int|\chi(\mathbf r)|^2e^{-i\mathbf p\cdot\mathbf r}d^3r,\qquad
n=\langle\chi|\chi\rangle .
$$

**Check:** $\tilde F(\mathbf 0)=0$, i.e. $\langle\tilde\phi_b|\chi\rangle=0$.

Everything else in the amplitude is unchanged: the projectile factors QT1, QC2, the
TSCWN terms and $C_{ab}$. The Gamow factor $|N_b|^2$ multiplies $\phi_b$ and $\tilde\phi_b$
in the same way.

## 3. Molecular orbital and orientation average: channel-by-channel orthogonalization

The code uses the orientation average (report, Eqs. 18–19). The OCE orbital is split
into parts $g$:

$$
\Phi=\sum_g R_g(r)\,S_{l_g m_g}(\hat r),\qquad
\langle|T|^2\rangle=\sum_g\frac{1}{2l_g+1}\sum_{\mu=-l_g}^{l_g}\big|T[\chi_{g\mu}]\big|^2,
\qquad \chi_{g\mu}=R_g\,Y_{l_g\mu}.
$$

This formula is exact because $T$ is linear in the orbital. Orthogonalizing to the full
$\Phi$ at each orientation would make $T$ non-linear in $\Phi$ (through $S$), so the
average would no longer reduce to this sum. That would require a numerical average over
the Euler angles.

We therefore orthogonalize **each channel to itself**:
$\chi\to\chi_{g\mu}$, $n\to n_g=\int R_g^2r^2dr$. This is exactly the atomic
prescription applied to every partial wave of the one-centre expansion, and it is
consistent with the incoherent channel sum already used. For CH$_4$ 1t$_2$ (pure p
part) it is identical to orthogonalizing an atomic 2p-like orbital.

## 4. Formulas used in the code

**Overlap** (once per angle and channel): $S_{g\mu}=(2\pi)^{3/2}F_{g\mu}(\mathbf 0)$.
- Analytic part: TFGA with $\mathbf{ad}=0$.
- TECWN part with $e^{-i\mathbf p\cdot\mathbf r}=1$: with the mirror-symmetric half
  grid, $E^+=2$ and $E^-=0$, so
  - $S_0=2\sum c_0$,
  - $S_\mu=\sum c^{(+\mu)}$,
  - $S_{-\mu}=(-1)^\mu\sum c^{(+\mu)}$.

  (Array `sov(icm,ia)`.)

**Density form factor.** With $e^{-i\mathbf p\cdot\mathbf r}=4\pi\sum_L(-i)^Lj_L(pr)\sum_M Y_{LM}(\hat p)Y^*_{LM}(\hat r)$,
and since $|Y_{l\mu}|^2$ contains only $Y_{L0}$ with $L$ even and $L\le 2l$:

$$
\frac{G_{g\mu}(\mathbf p)}{n_g}=\sum_{L=0,2,4}\underbrace{\frac{4\pi}{(2\pi)^{3/2}}\frac{(-i)^L}{n_g}A_{l\mu L}\sqrt{\tfrac{2L+1}{4\pi}}}_{\texttt{cgo(L/2,icm)}}\;
\underbrace{\int_0^\infty R_g^2(r)\,j_L(pr)\,r^2dr}_{\texttt{aiq(L/2,g)}}\;P_L(\cos\theta_p),
$$

$$
A_{l\mu L}=\int|Y_{l\mu}(\hat r)|^2Y_{L0}(\hat r)\,d\Omega .
$$

- $G$ is real and depends only on $|\mathbf p|$ and $\theta_p$.
- Radial integral: Gauss–Legendre on $[0,30]$ a.u., 15 intervals $\times$ `nleg` points.
  It is computed once per parallel task, since $|\mathbf p|$ is fixed per task.
- $A_{l\mu L}$: Gauss–Legendre in $\cos\theta$ (exact for these polynomials).
- $j_0,j_2,j_4$: closed forms for $x\ge1$, power series for $x<1$ (`SPHJ`).

**Amplitude:**
$$
T_{g\mu}=\int d^3p\;\big[F_{g\mu}(\mathbf p)-S_{g\mu}G_{g\mu}(\mathbf p)/n_g\big]
\big[Q_{C1}(Q_{T1}+Q_{W1})-Q_{T2}(Q_{C2}+Q_{W2})\big].
$$

The TDCS is unchanged:

$$
\sigma^{(3)}=2\,a_{Ns}\,N_aN_bN_{ab}\sum_{g\mu}\frac{|T_{g\mu}|^2}{2l_g+1}.
$$

## 5. Validation

- `iorth=0` reproduces the previous code exactly (same operations).
- $G(\mathbf p)/n$ from the formula above vs brute-force 3D integration of
  $(2\pi)^{-3/2}n^{-1}\int|R\,Y_{l\mu}|^2e^{-i\mathbf p\cdot\mathbf r}d^3r$ at
  $\mathbf p=(0.7,-0.4,0.9)$, all channels: relative difference $<6\times10^{-8}$
  (printed by the program at start).
- By construction $\tilde F(\mathbf 0)=0$.

**Overlaps $|\langle\phi_b|\chi_{g\mu}\rangle|$** at $\theta_e=0$, $E_e=74$ eV (printed by the program):

| Orbital | Overlaps |
|---|---|
| 2a$_1$ | s: 0.37; p, d: $\le$ 0.02 |
| 3a$_1$ | s: 0.10; p: 0.18–0.20; d: $\le$ 0.01 |
| 1e | d: 0.003–0.03; p: 0.028 |

The largest overlaps are in the channels where the short-range term overshot.

## 6. Limits of the method

- Orthogonalization removes the spurious part of $T^0_{dir}$. It does **not** correct
  the short-range phase of the variable-charge Coulomb wave. That needs a distorted wave:
  a numerical solution of $[-\tfrac12\nabla^2-Z_{ion}(r)/r]\phi=E\phi$, the next step.
- Channel-by-channel orthogonalization is an approximation for a molecule (§3). The exact
  alternative is to orthogonalize to the full orbital at each orientation and average
  numerically over the Euler angles.
- Exchange ($T_{exc}$) is still neglected; it is small for $E_e=74$ eV, $E_s=500$ eV.
