# Independent GetDP complete-electrode verification

Implemented and run **26 October 2026**. [Code](../examples/getdp_cem.py) · [native formulation](../examples/getdp_cem.pro) · [nine-case release](../data/synthetic_getdp/README.md) · [red–green record](tdd-getdp.md).

GetDP now supplies a second FEM implementation of the complete-electrode model (CEM). It assembles and solves its own system and evaluates its own fields, electrode currents and power integrals. The standalone adapter imports NumPy and Gmsh, **not DOLFINx**. The comparison driver imports both backends, shares meshes and SI inputs, and samples the DOLFINx fields at GetDP's tetrahedral barycenters. This checks implementation independence on matched discretizations; it does not establish anatomical or experimental accuracy.

## Formulation and conventions

Tissue potential uses GetDP's `BF_Node` basis, supplemented by hierarchical `BF_Node_2E` functions for P2. Each metal electrode has its own region-constant `BF_Region` voltage. The contact law permits tissue potential to vary beneath an equipotential metal electrode. GetDP's [function-space and integration documentation](https://getdp.info/doc/texinfo/getdp.html) describes these primitives.

The independently written weak form is

\[
\int_\Omega \sigma\nabla u\cdot\nabla v\,dx
+\sum_l\int_{e_l}\frac{(u-U_l)(v-V_l)}{z_l}\,ds
=\sum_l I_l V_l.
\]

Currents are inward-positive, sum to zero, and contact impedances have units Ω·m². One electrode voltage fixes the gauge. Scalar or piecewise symmetric positive-definite tensors are supported. Unlisted exterior facets are insulated. Three floating-electrode patterns exercise reciprocity, voltage-gauge changes, polarity and zero current. In that test, tag 199 deliberately becomes a third metal electrode; in the released two-electrode benchmarks, it remains insulation.

The mesh geometry is linear tetrahedral, with tagged triangular boundaries. Four-point tetrahedral and six-point triangular Gaussian rules integrate the P2 bulk and contact terms for the piecewise-constant materials. The native solve uses PETSc LU through GetDP; our DOLFINx reference uses SciPy's serial direct solve. No matrix or potential vector is transferred between solvers.

## Evidence and its scope

| Check | Result in the checked-in macOS release |
|---|---|
| Full-face homogeneous / layered analytical voltages | 51.4583333 / 79.5833333 mV at 1 mA |
| Nine same-mesh cases: analytical slabs, finite patches, off-diagonal tensor, adaptive mesh | Maximum relative voltage difference **1.20×10⁻¹²** |
| Current field at common cell barycenters | Maximum relative sample RMS difference **1.19×10⁻¹²** |
| Integrated electrode currents | Maximum prescribed-current discrepancy **4.69×10⁻¹⁷ A** |
| Adaptive mesh, independent GetDP P2 solve | **73.9996092 Ω** |
| Primal/complementary resistance interval, recomputed after XDMF reload | **73.9996092–74.6038179 Ω** |
| Independently generated finest global Gmsh mesh, GetDP P2 | **74.0073682 Ω**, inside that interval |

The three globally generated meshes give 73.4760846, 73.7505454 and 74.0073682 Ω. Coarser Galerkin results fall below the adaptive lower bound, as expected for less accurate primal approximations; they are not required to lie inside a finer-mesh interval. The GetDP solve on the *same adaptive mesh and order* reproduces its primal endpoint. The upper endpoint comes from DOLFINx's admissible RT reconstruction, not from an independent GetDP dual solve.

Sample RMS means `norm(GetDP − DOLFINx) / norm(DOLFINx)` over unweighted tetrahedral barycenters. It is **not** a volume-weighted L2 norm, a pointwise bound, or an estimate of error against the exact solution. Agreement near roundoff can coexist with the released 0.8165% complementary-energy gap. The resistance interval is numerical, subject to the assumptions and floating-point qualifications in the [conservative-current model](fem-conservative-current.md); it does not bound MRI phase.

The XDMF→MSH exporter maps topology entities to geometric DOFs explicitly, retains coordinates at round-trip double precision and orients tetrahedra positively. Tests compare every cell/facet's geometric connectivity and physical label, check both electrode areas, and re-solve to verify voltage and current signs. The output can be used by GetDP without XDMF or DOLFINx.

## Reproduce

Install the [isolated FEM environment](fem-complete-electrode.md#reproduce) for comparisons. The standalone `solve_getdp` adapter needs only NumPy, Gmsh and a GetDP executable. Set `GETDP_EXECUTABLE` to the executable path, or put `getdp` on `PATH`.

Linux CI downloads the official **3.5.0** archive from [GetDP](https://getdp.info/bin/Linux/getdp-3.5.0-Linux64c.tgz), verifies SHA-256 `d3c28fa18f20d6147b4c7367d4dd802e9f7ddb58c608688bbb71919dbca8041d`, and checks `getdp -info` before running tests. The exact installation steps are in [the workflow](../.github/workflows/benchmarks.yml).

The checked-in numerical release was generated on Apple silicon using official **4.0.0-git-259f4c38**, built 12 October 2026, PETSc 3.17.4. Its [macOS ARM snapshot archive](https://getdp.info/bin/macOS/getdp-git-MacOSARMc.tgz) had SHA-256 `31b106b32ea88632764f681310245f3e456bdd521aa781c51f861bd1b8e0019c`. That URL is mutable: the recorded hash identifies the downloaded archive, but does not guarantee that a later download retrieves it. Each release additionally records the executable hash and full native `-info` output. No third-party binary is checked in.

```sh
python -m unittest discover -s tests -p 'test_getdp*.py' -v
python examples/getdp_benchmark.py --output /tmp/synthetic-getdp --plot
```

For a native replay, enter a case directory such as `adaptive_p2` and run:

```sh
getdp problem.pro -msh ../meshes/adaptive.msh -solve CEM -pos Export -ksp_type preonly -pc_type lu
```

Choose the case's mesh listed in `parameters.json`; for example both tensor cases use `finite_patch_0.msh`. Native result files, numerical tables, sampled comparison arrays, logs, source hashes and input hashes are retained. Logs include elapsed time and timestamps, so byte-identical regeneration is not promised. The tests require numerical agreement without adjusting tolerances between macOS development GetDP and Linux stable GetDP.

The new nine tests bring the suite to **71 distinct tests**: 37 spectroscopy/FV, 25 existing DOLFINx FEM and nine GetDP/comparison tests. [Linux run 40455074002](https://github.com/m9h/neuroelectric-mri/actions/runs/40455074002) passed both jobs, including the independent solver. Complete-circuit B, observable-level error control, anatomical geometries, MPI scaling and validated derivatives remain separate gates in the [roadmap](fem-roadmap.md).
