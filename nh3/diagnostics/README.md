# Diagnostic: is the 3CWZ ejected wave a good solution of Eq. (11)?

- `drv.f`: prints Z_ion(r) from ZSCR (link with the subroutines of `nh3_opt_claude.for`
  from `SUBROUTINE TARGET` to the end): `gfortran drv.f lib.f`, output in fort.11-13.
- `pw.py`: Python (numpy, scipy, mpmath). It solves the radial equation
  u'' + [k^2 - l(l+1)/r^2 + 2 Z_ion(r)/r] u = 0 at E_e = 74 eV.
  - It extracts the short-range phase shifts δ_l (beyond Coulomb Z=1).
  - It compares the partial-wave content of the exact wave, the variable-charge 3CWZ
    wave and the Z=1 Coulomb wave for r < 2.5 a.u.
  - It gives the elastic e⁻ + NH₃⁺ cross section relative to Rutherford (Z=1).
The results are in `../theory_distorted_wave.md`.
