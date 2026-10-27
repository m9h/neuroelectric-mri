# 🟣 Independent GetDP × DOLFINx CEM benchmark

**Synthetic only.** Nine independently assembled comparisons cover analytical homogeneous/layered slabs, three globally refined finite-electrode meshes, P1/P2, an off-diagonal SPD tensor and the exported adaptive mesh. [Numerics and reproduction](../../docs/getdp-verification.md) · [red–green record](../../docs/tdd-getdp.md).

![Independent solver comparison and resistance interval](summary.png)

GetDP and DOLFINx differ by at most **1.20×10⁻¹²** in relative electrode voltage and **1.19×10⁻¹²** in cell-center current sample RMS in this release. These measure agreement between implementations, not error against the continuum solution. The adaptive GetDP resistance is **73.9996092 Ω**; the recomputed primal/complementary interval is **73.9996092–74.6038179 Ω**. The independently generated finest global mesh gives **74.0073682 Ω**. The numerical interval is not an MRI-phase or pointwise-current bound.

| File | Content |
|---|---|
| `meshes/` | Original tagged Gmsh meshes plus the tested XDMF→MSH adaptive conversion |
| `metrics.csv` | Voltage, separate bulk/contact powers, current residuals, sampled field differences and adaptive resistance interval |
| `parameters.json` | Exact conductivity/current/contact/gauge inputs, mesh used by each case, SI and sampling conventions |
| `<case>_p<degree>/problem.pro`, `parameters.pro` | Independent native GetDP formulation and case inputs |
| `<case>_p<degree>/problem.pre`, `problem.res` | Native preprocessed problem and solution coefficients; replay with the specified mesh |
| `<case>_p<degree>/*.txt` | Native field/power/current/voltage tables and complete solver build information |
| `<case>_p<degree>/comparison.npz` | Common coordinates and each solver's potential, E and J; read with `allow_pickle=False` |
| `<case>_p<degree>/solver.log` | Native solver log |
| `manifest.json` | SHA-256 hashes of input meshes, generated artifacts, scientific source files and GetDP executable, plus Python/package/platform versions |

`comparison.npz` fields use m, V, V/m and A/m². The reported RMS is an **unweighted sample norm at tetrahedral barycenters**, not a volume L2 norm. GetDP current is the primal `−σ∇u`; it is not the equilibrated RT current. The dual interval is recomputed with the existing DOLFINx RT implementation. Small platform/version floating-point differences and different log timestamps are expected.

Generated on Apple silicon with GetDP **4.0.0-git-259f4c38** and the pinned DOLFINx development environment; Linux CI independently runs the same assertions with GetDP **3.5.0**. [Archive hashes and installation details](../../docs/getdp-verification.md#reproduce). No solver binary, B field, scanner data or experimental validation is included. Generated numerical artifacts use **CC0-1.0**; third-party software retains its own license.
