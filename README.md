# PYTHIA low-photon-energy cutoff study

This package documents the low-energy cutoff observed in a PYTHIA 8.3
photoproduction study intended as a noncollective baseline for inclusive
photonuclear collisions at RHIC.

## Physics setup

- Beam configuration: 10 GeV electron on 100 GeV proton.
- Photons from the electron are generated with the equivalent-photon
  approximation using `PDF:lepton2gamma = on`.
- DIS is disabled with `WeakBosonExchange:ff2ff(t:gmZ) = off`.
- The October 8 generator configuration uses `Photon:Q2max = 0.01`,
  `Photon:Wmin = 5.0`, `Photon:Wmax = -1.0`, and
  `Photon:ProcessType = 1`.
- `SoftQCD:nonDiffractive` is enabled, while the inclusive `SoftQCD`,
  `HardQCD`, and `PhotonParton` switches are first disabled explicitly.
- Multiparton interactions, initial- and final-state radiation,
  hadronization, and colour reconnection are enabled in this configuration.
- The analysis reconstructs the exchanged photon from the electron beam line
  and selects `0 <= Q2 < 0.01 GeV^2` and `0 < E_gamma < 5 GeV`.
- The generated photon spectrum can subsequently be reweighted toward the
  neutron-tagged photonuclear flux stored in `inputs/noon_flux.root`. Event
  weights are intentionally set to one in the supplied cutoff-test analysis.

For a head-on photon-proton collision at high energy,

```text
W_gamma-p^2 approximately equals 4 E_gamma E_p.
```

With `E_p = 100 GeV` and `W_gamma-p = 5 GeV`, the corresponding photon
energy is approximately 62.5 MeV. This agrees with the observed low-energy
edge and is correlated with the cutoff of the raw final-state multiplicity
near four particles.

## Package contents

- `code/simrun15mb_ep.C`: October 8 STAR generator configuration. The relevant beam and
  photoproduction settings are around lines 270-330, with the principal
  PYTHIA switches on lines 309-330.
- `code/analysis.cxx`: analysis executable that fills `egamma`,
  `h_npart_raw`, and `h_npart_tof_matched`. The exchanged-photon
  reconstruction and event selection are around lines 323-373.
- `code/reweight_egamma.C`: construction of photon-energy weights from the
  nOOn flux.
- `code/make_diagnostic_plots.C`: creates a two-panel photon-energy and raw
  multiplicity diagnostic figure from the representative ROOT output.
- `inputs/noon_flux.root`: neutron-tagged photonuclear flux input.
- `inputs/reweighting_factors_ep2x100_test2.txt`: representative weight table (obtained using n00n model).
- `inputs/egamma_less1p1gev_ep10x100_test17-18_eta2_01.root`: ROOT file used as the input to the diagnostic
  plotting macro.
- `outputs/Egamma_reweighted.pdf`: existing photon-spectrum reweighting plot.

## Reproduce the diagnostic figure

From the package directory, run:

```bash
root -l -b -q 'code/make_diagnostic_plots.C()'
```

This creates `outputs/low_energy_cutoff_diagnostics.pdf` and
`outputs/low_energy_cutoff_diagnostics.png`.

By default, the macro reads
`inputs/egamma_less1p1gev_ep10x100_test17-18_eta2_01.root`
