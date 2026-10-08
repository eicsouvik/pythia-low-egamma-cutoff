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

- `code/simrun15mb_ep.C`: October 8 STAR generator configuration copied from
  `Code/upc2pc/model/macros/10_8_2026/simrun15mb_ep.C`. The relevant beam and
  photoproduction settings are around lines 270-330, with the principal
  PYTHIA switches on lines 309-330.
- `code/analysis.cxx`: analysis executable that fills `egamma`,
  `h_npart_raw`, and `h_npart_tof_matched`. It is copied from
  `Code/upc2pc/model/macros/10_8_2026/analysis.cxx`; the exchanged-photon
  reconstruction and event selection are around lines 323-373.
- `code/reweight_egamma.C`: construction of photon-energy weights from the
  nOOn flux.
- `code/make_diagnostic_plots.C`: creates a two-panel photon-energy and raw
  multiplicity diagnostic figure from the representative ROOT output.
- `inputs/noon_flux.root`: neutron-tagged photonuclear flux input.
- `inputs/reweighting_factors_ep2x100_test2.txt`: representative weight table.
- `inputs/egamma_less1p1gev_ep10x100_test17-18_eta2_01.root`: an earlier
  representative analyzed ROOT file used as the input to the diagnostic
  plotting macro. It was not regenerated with the October 8 code; replace it
  with a new output before making quantitative comparisons.
- `outputs/Egamma_reweighted.pdf`: existing photon-spectrum reweighting plot.

## Reproduce the diagnostic figure

From the package directory, run:

```bash
root -l -b -q 'code/make_diagnostic_plots.C()'
```

This creates `outputs/low_energy_cutoff_diagnostics.pdf` and
`outputs/low_energy_cutoff_diagnostics.png`.

By default, the macro reads
`inputs/egamma_less1p1gev_ep10x100_test17-18_eta2_01.root`. A different ROOT
file can be supplied as the first macro argument.

## Questions for the PYTHIA authors

1. Is the approximately 5 GeV lower bound on the gamma-p invariant mass a
   hard limit of the lepton-to-photon flux or photon-hadron initialization in
   this PYTHIA configuration?
2. Can lower-energy photons be generated consistently within the same EPA
   setup with `SoftQCD:nonDiffractive`, or should a different
   soft-photoproduction configuration be used below this limit?
3. Are additional `PhaseSpace`, photon PDF, VMD/GVMD, direct/resolved, or
   minimum-mass settings required?
4. Is a hybrid sample needed to cover low-mass hadronic final states that are
   outside the perturbative photon-parton configuration?

## October 8 analysis updates

The supplied analysis no longer assumes a fixed event-record index for the
exchanged photon. It identifies a unique non-final-state photon connected to
the incoming electron beam line, computes `Q2 = -q.M2()`, and fills the photon
spectrum using the standard weighted ROOT call
`egamma.Fill(E_gamma * 1e3, event_weight)`.
# pythia-low-egamma-cutoff
