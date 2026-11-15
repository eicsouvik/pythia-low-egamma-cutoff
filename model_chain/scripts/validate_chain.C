#include "TFile.h"
#include "TH1.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

template <class T>
T *require(TFile &file, const char *name)
{
  T *object = dynamic_cast<T *>(file.Get(name));
  if (!object) throw std::runtime_error(TString::Format(
      "Missing '%s' in %s", name, file.GetName()).Data());
  return object;
}

void requireBranch(TTree *tree, const char *name)
{
  if (!tree->GetBranch(name)) throw std::runtime_error(TString::Format(
      "Missing branch '%s' in tree '%s'", name, tree->GetName()).Data());
}

} // namespace

void validate_chain(const char *pythiaPath, const char *neutronPath)
{
  try {
    std::unique_ptr<TFile> pythia(TFile::Open(pythiaPath, "READ"));
    if (!pythia || pythia->IsZombie())
      throw std::runtime_error("Cannot open weighted PYTHIA output");
    TNamed *flux = require<TNamed>(*pythia, "noonFluxWeightApplied");
    if (!TString(flux->GetTitle()).BeginsWith("yes"))
      throw std::runtime_error("PYTHIA output does not record applied n00n weights");
    require<TNamed>(*pythia, "noonDetectorResponse");
    TTree *interfaceTree = require<TTree>(*pythia, "noonEvents");
    for (const char *branch : {"eventIndex", "egammaGeV", "nTofMatched",
                               "multClass", "generatorWeight",
                               "noonFluxWeight", "totalWeight"})
      requireBranch(interfaceTree, branch);

    std::unique_ptr<TFile> neutron(TFile::Open(neutronPath, "READ"));
    if (!neutron || neutron->IsZombie())
      throw std::runtime_error("Cannot open n00n neutron output");
    TParameter<int> *complete = require<TParameter<int>>(*neutron, "runComplete");
    if (complete->GetVal() != 1)
      throw std::runtime_error("n00n output is an incomplete checkpoint");
    TNamed *neutronFlux = require<TNamed>(*neutron, "noonFluxWeightApplied");
    if (!TString(neutronFlux->GetTitle()).BeginsWith("yes"))
      throw std::runtime_error("n00n output did not preserve applied flux metadata");
    TNamed *neutronInput = require<TNamed>(*neutron, "inputFile");
    if (TString(neutronInput->GetTitle()) != TString(pythiaPath))
      throw std::runtime_error("n00n input provenance does not match the PYTHIA file being validated");
    require<TTree>(*neutron, "neutronEvents");

    for (int activity = 0; activity < 2; ++activity) {
      TH1 *before = require<TH1>(*neutron,
          Form("hEgammaInput_c%d", activity));
      TH1 *after = require<TH1>(*neutron,
          Form("hEgammaAfterburner_c%d", activity));
      require<TH1>(*neutron, Form("hNeutronTotal_c%d", activity));
      if (before->GetNbinsX() != after->GetNbinsX())
        throw std::runtime_error("Photon-energy closure binning mismatch");
      for (int bin = 0; bin <= before->GetNbinsX() + 1; ++bin) {
        const double difference = std::abs(
            before->GetBinContent(bin) - after->GetBinContent(bin));
        const double scale = 1.0 + std::abs(before->GetBinContent(bin));
        if (difference > 1.0e-12 * scale)
          throw std::runtime_error(TString::Format(
              "Photon-energy closure failed for class %d, bin %d",
              activity, bin).Data());
      }
    }

    std::cout << "VALIDATION PASSED\n"
              << "  n00n flux weighting metadata: yes\n"
              << "  STAR-like detector response metadata: present\n"
              << "  interface branches: present\n"
              << "  n00n runComplete: 1\n"
              << "  PYTHIA-to-n00n photon-energy closure: passed"
              << std::endl;
  } catch (const std::exception &error) {
    std::cerr << "VALIDATION FAILED: " << error.what() << std::endl;
    gSystem->Exit(2);
  }
}

