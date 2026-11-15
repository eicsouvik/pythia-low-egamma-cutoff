# Provenance and scope

The `noon_original/NeutronGenerator.*` implementation and nuclear-data tables are
derived from the public [mbroz84/noon](https://github.com/mbroz84/noon)
repository at commit `2b5056c`. The interface and QA macros in this package
are analysis-specific additions. This package does not relicense the upstream
code; users should consult the upstream project and authors before
redistribution or publication.

The PYTHIA generator and analysis code records the configuration used in this
study. Running it requires the STAR software environment and the corresponding
event-tree interface (`genevents`). The generated event files are not included
because they are large and installation-specific.

The provided ROOT input `n00n_flux_from_Ashik_da/noon_flux.root` and weight
table are research
inputs for reproducibility of this analysis workflow. They are not a general
purpose photonuclear event sample.
