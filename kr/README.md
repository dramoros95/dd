# Kr 4p (e,3e) double ionization: optimized BBK code with variable charge

`kr3e_opt.for` is an optimized and corrected version of `original/DI-Kr-4p-BBK-zv-3P.for`.
- **Physics:** unchanged. Fast projectile in first Born; BBK final state for the two slow
  electrons, with variable charge Z(r); C12 handled by Fourier transform over p; Kr²⁺ 4p⁴ ³P.
- **Implementation:** reorganized for speed and accuracy.
- **Wave functions:** the Clementi tables (4.1–4.3).

## Compile and run
```bash
gfortran -O3 -march=native -fopenmp -ffixed-line-length-132 kr3e_opt.for -o kr3e_opt
# or: ifx -O3 -xHOST -qopenmp -extend-source 132 kr3e_opt.for -o kr3e_opt
export OMP_NUM_THREADS=4
echo "1 73" | ./kr3e_opt        # theta_1 = 5*(index-1) deg, output kr3e_3P.dat
```

Switches at the top of the main program:

| Switch | Value | Meaning |
|---|---|---|
| `iwf` | 1 | Clementi tables |
| | 2 | Data of the original code |
| `izv` | 1 | Variable charge Z(r) |
| | 0 | Constant Z = 2 |
| `lmax`, `dz`, `ndr`×`nlr` | | Partial waves and radial grid of the Z(r) correction |

The kinematics (Eₛ, E₁, E₂, θₛ, θ₂, step in θ₁) and the p grid are set just below the switches.

## What was changed, and why
1. **Z(r) term integrated with one point.** The original `tdN` used `nlg=nlg1=1`: one point
   in r, θ and φ. The term is now computed by partial waves (`PWDL`, `JTAB`, `CCOEF`, `PWQ3`):
   - the radial integrals are tabulated once;
   - the angular integration is exact (Gaunt coefficients).
2. **Typo in the 4p coefficient.** `Fr4p`, `df4p` and `qkr4p` used 0.003172 instead of
   0.03172, which gave a 4p norm of 1.0041. The 4p data now come from one place.
3. **Missing cross terms in Z(r).** `az2p`, `az3p` and `az4p` dropped the 2p–4p cross terms
   (summed into `av2s4s`). The new general `ZKR` routine handles any Slater basis.
4. **Clementi tables.** Two problems in table 4.1:
   - 3S column, 9th row: 0.02180 should read 0.00218. With 0.02180 the 3s norm is 0.990
     and ⟨3s|4s⟩ = 0.019.
   - The 2S column is not orthonormal (norm 0.982, ⟨1s|2s⟩ = 0.019), probably several
     misprints. It enters only Z(r), through its density, and is renormalized. **Check it
     against the source.**

   All orbitals are normalized at run time and their norms are printed.
5. **Speed:**
   - `TCP3` computes m = −1, 0, 1 in one call;
   - electron 2 (fixed k₂) is computed once for all angles;
   - OpenMP over the p grid;
   - the 960 MB of static arrays are removed.

## Validation (`tests/`)
| Test | Result |
|---|---|
| `t_tcp3_ylm.f`: TCP3 vs the original `tcp` | identical (0 difference) |
| `t_tcp3_ylm.f`: YALL vs `ylm` | 1e-10 (π truncated in `ylm`) |
| `t_pwcor.f`: partial-wave Z(r) term vs brute-force 3D integral | ≤ 6e-8 relative |
| `t_regression.py`: optimized (`iwf=2`, `izv=0`) vs original with the Z(r) term removed | ratio 1.00000248 at every angle (= 1/norm² of the original 4p data) |

The analytic Z = 2 part (`tcp`) and the numerical Z(r) part (`tdN`) of the original code use
the same conventions: they agree to 6 digits for Z = 2.

## Cost
| Run | Time |
|---|---|
| Original (Z(r) term with 1 point) | ≈ 45 s per angle on 1 core (≈ 55 min for 73 angles) |
| Optimized (converged Z(r) term) | ≈ 5 s per angle on 4 cores (≈ 6 min for 73 angles) |
