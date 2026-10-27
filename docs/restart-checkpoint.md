# Restart checkpoint — 26 October 2026

Saved at the user's request before a system reboot.

## Completed

- Independent native GetDP CEM, P1/P2, scalar/SPD tensors, contact currents/power, floating-electrode reciprocity, gauge and polarity checks.
- Nine-case GetDP/DOLFINx comparison and tagged adaptive XDMF→MSH export with geometric round-trip/current-sign tests.
- Synthetic release in `data/synthetic_getdp/`, including native solver inputs/results, comparison arrays, plots and SHA-256 provenance.
- Red–green commits: `83f0e83` → `0ff5cce`; `41d879f` → `1a68f8d`.
- Local checks: 37 spectroscopy/FV + 25 DOLFINx + nine GetDP/comparison tests pass (71 distinct methods).
- Both GitHub Actions jobs passed at `1a68f8d`: [run 40455074002](https://github.com/m9h/neuroelectric-mri/actions/runs/40455074002). Linux uses checksum-pinned GetDP 3.5.0; the local release uses 4.0.0-git-259f4c38. Documentation/data changes after that run have not changed solver code.
- Full interpretation and installation details: [GetDP verification](getdp-verification.md). The README and roadmaps now identify complete-circuit B as the next ticket.

## Resume

1. Check `git status`, remote main and the latest Actions run for this documentation/data checkpoint.
2. Continue the [complete-circuit B/phase ticket](fem-roadmap.md#next-ticket-complete-circuit-magnetic-field-and-mr-phase), keeping the user's red–green TDD requirement. Start with tests and a fully specified closed current path; electrode connections and the return lead must be included.
3. Retain the isolated FEM environment. The local ignored executable is `.tools/getdp-git-MacOSARM/bin/getdp`; `.venv` contains the FEM packages. Native local commands use `env OMPI_MCA_btl=self XDG_CACHE_HOME=/private/tmp/neuroelectric-fem-cache .venv/bin/python ...`. Temporary caches/logs may disappear after reboot; all scientific outputs and provenance needed for this milestone are in the repository.

Scope still open: complete-circuit B, near-field quadrature, MR observable error control, broader geometry/tensor studies, PETSc/MPI, local RT reconstruction, MFEM throughput evaluation, verified derivatives and measured Prisma phantom data. No hardware or scanner was operated, and no external message was sent. The existing mrs-jax work remains on its separate PR; this checkpoint does not change it.
