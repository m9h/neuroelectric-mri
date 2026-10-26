# Complete-electrode FEM reference

Implemented in [fem_tank.py](../examples/fem_tank.py), with [synthetic data](../data/synthetic_fem/README.md) and a [red–green record](tdd-fem.md). This is a DOLFINx tetrahedral FEM implementation independent of the earlier Cartesian finite-volume tank.

## Physical convention and weak form

Use metres, volts, amperes and S/m. Let u be the tissue potential, Uₑ the equipotential metal electrode voltage, σ a symmetric positive-definite conductivity tensor, and zₑ the strictly positive area-specific contact impedance in **Ω·m²**. Currents Iₑ are positive **into** the tank and sum to zero. This differs from the older FV boundary ledger, which is outward-positive.

In the bulk, J = −σ∇u and ∇·J = 0. At electrode e, the inward current density is (Uₑ−u)/zₑ. Unlisted boundary tags are insulated; a listed electrode with Iₑ=0 is floating and can carry local current of both signs.

The complete electrode model has weak form

\[
\int_\Omega \nabla v\cdot\sigma\nabla u\,dV
+\sum_e z_e^{-1}\int_{E_e}(u-U_e)(v-V_e)\,dS
=\sum_e I_e V_e.
\]

This is the classical CEM; see [Somersalo, Cheney & Isaacson (1992)](https://doi.org/10.1137/0152060). We assemble the field/electrode block matrix with DOLFINx/UFL and fix one electrode voltage to zero. Its redundant current equation is retained in the post-solve residual check. The ideal z=0 limit needs a constrained formulation and is rejected rather than approximated silently.

CG1/CG2 represents u; piecewise-constant scalar/tensor σ is assigned by material tags. SciPy solves the small system in serial. The returned E and J are DG0/DG1 gradient reconstructions. **They are not locally conservative H(div) fields.** Integrated electrode currents use the CEM contact law, not an uncorrected gradient trace. The independent energy check is

\[
\sum_e U_e I_e
=\int_\Omega \nabla u\cdot\sigma\nabla u\,dV
+\sum_e z_e^{-1}\int_{E_e}(u-U_e)^2\,dS.
\]

## Geometry and validation

OpenCASCADE Boolean fragments conform the mesh to material interfaces and electrode edges. Tags 1/2 are materials, 101/102 are electrodes on x−/x+, and 199 is the remaining exterior surface. The default tank is 120×80×60 mm. A patch fraction of 0.5 produces 40×30 mm electrodes, one quarter of the full-face area. Finite edges have mesh size h/8 grading to h over distance h; this is prescribed refinement.

The eleven FEM tests cover CAD volume/area and MSH round-trip, an equal-area electrode-tag regression, closed-form contact drops and power, layered series resistance and independent FV comparison, a rotated off-diagonal conductivity tensor, gauge/polarity/off, floating electrodes and reciprocity, invalid inputs, h/p voltage refinement, exported fields, hashes and repeat generation.

At 1 mA with full-face contacts z=(0.002, 0.005) Ω·m², the homogeneous 0.5 S/m voltage is **51.458333 mV** and the equal-thickness 0.2/0.8 S/m layered voltage is **79.583333 mV**. Both polynomial orders agree with these references to numerical precision. The finite-patch CG2 voltages are **73.4761, 73.7505 and 74.0074 mV** across the published mesh levels. The last change is **0.347%**. This is a successive-resolution observation, not a bound on the remaining voltage or local-current error. The finest CG1 voltage still differs from CG2 by about **3.55%**.

## Reproduce

Use a separate Python 3.12 environment; the pinned FEniCS development wheels support macOS 15+ and Linux with glibc 2.34+. Linux also needs `libGLU` for the Gmsh wheel. See [official installation instructions](https://github.com/FEniCS/dolfinx/blob/main/README.md).

```sh
python3.12 -m venv .venv-fem
.venv-fem/bin/python -m pip install -r requirements-fem.txt
.venv-fem/bin/python -m pip install matplotlib
.venv-fem/bin/python -m unittest discover -s tests -p 'test_fem*.py' -v
.venv-fem/bin/python examples/fem_benchmark.py --output /tmp/synthetic-fem --plot
```

The FEM CI job checks imports explicitly, runs all eleven tests without skips and exports an artifact. The spectroscopy/FV job runs the previous 37 tests in its own environment; optional FEM tests skip there. These are 48 distinct tests across the two environments.

Remaining work is specified in the [FEM roadmap](fem-roadmap.md): conservative currents, error estimation, independent GetDP CEM, complete-circuit magnetics, PETSc/MPI, derivatives and measurements. There is no conductivity inversion, anatomical validation or scanner acquisition in this release.
