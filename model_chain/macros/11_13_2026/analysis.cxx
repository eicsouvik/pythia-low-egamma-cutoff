#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
//#include <stdio.h>
#include <stdlib.h>
//#include "iostream.h"
//#include "iomanip.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"

#include "TChain.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TParticlePDG.h"
#include "TNamed.h"
#include "TLeaf.h"
#include "TKey.h"
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "genevents.h"
using namespace std;

//define functions
bool isPrimaryParticle(const int id);
bool isChargedParticle(const int id);
void findMotherParticle(genevents *tree);

// Mother indices in the STAR primary record already include generator offsets.
// Follow only same-sign electron copies: photons from hadron decays or other
// final-state radiation must never stand in for the exchanged beam photon.
bool electronBeamLine(const genevents *t, int index, int id)
{
  for (int depth = 0; depth < t->mNumParticles; ++depth) {
    if (index <= 0 || index >= t->mNumParticles ||
        t->mParticles_mId[index] != id) return false;
    if (t->mParticles_mStatus[index] == 4) return true;
    int next = -1;
    for (int m = 0; m < 2; ++m) {
      int parent = t->mParticles_mMother[index][m];
      if (parent > 0 && parent < t->mNumParticles && parent != index &&
          t->mParticles_mId[parent] == id) {
        if (next != -1 && next != parent) return false;
        next = parent;
      }
    }
    index = next;
  }
  return false;
}

bool exchangedPhoton(const genevents *t, TLorentzVector &q)
{
  int photon = -1;
  for (int j = 0; j < t->mNumParticles; ++j) {
    if (t->mParticles_mId[j] != 22 || t->mParticles_mStatus[j] == 1) continue;
    bool fromBeam = false;
    for (int m = 0; m < 2; ++m) {
      int parent = t->mParticles_mMother[j][m];
      if (parent > 0 && parent < t->mNumParticles &&
          abs(t->mParticles_mId[parent]) == 11 &&
          electronBeamLine(t, parent, t->mParticles_mId[parent])) fromBeam = true;
    }
    if (!fromBeam) continue;
    if (photon != -1) return false; // Ambiguous records are not guessed.
    photon = j;
  }
  if (photon < 0) return false;
  q.SetPxPyPzE(t->mParticles_mPx[photon], t->mParticles_mPy[photon],
              t->mParticles_mPz[photon], t->mParticles_mEnergy[photon]);
  return true;
}


// A missing or ambiguous incoming beam proton is never guessed.
bool incomingProton(const genevents *t, TLorentzVector &p)
{
  int index = -1;
  for (int j = 0; j < t->mNumParticles; ++j) {
    if (t->mParticles_mId[j] != 2212 || t->mParticles_mStatus[j] != 4) continue;
    if (index != -1) return false;
    index = j;
  }
  if (index < 0) return false;
  p.SetPxPyPzE(t->mParticles_mPx[index], t->mParticles_mPy[index],
              t->mParticles_mPz[index], t->mParticles_mEnergy[index]);
  return TMath::Finite(p.Px()) && TMath::Finite(p.Py()) &&
         TMath::Finite(p.Pz()) && TMath::Finite(p.E()) && p.E() > 0;
}

// Optional, externally validated UPC weight and *physical flux density*.
// Both objects use x = laboratory E_gamma [MeV]; flux is dN/dE [1/MeV].
// A weight ratio is not a flux density and cannot determine missing flux.
// No extrapolation or clipping is allowed.  A validated exact-bin text table
// may be supplied as the fourth optional command-line argument.
struct NoonFluxInput {
  TFile *file;
  TObject *weights;
  TObject *density;
  NoonFluxInput() : file(0), weights(0), density(0) {}
  ~NoonFluxInput() { if (file) { file->Close(); delete file; } }
  static bool validate(TObject *o) {
    TH1 *h = dynamic_cast<TH1*>(o);
    TGraph *g = dynamic_cast<TGraph*>(o);
    if (h && h->GetDimension() == 1) {
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        if (!TMath::Finite(h->GetBinContent(b)) || h->GetBinContent(b) < 0 ||
            !TMath::Finite(h->GetXaxis()->GetBinLowEdge(b)) ||
            !TMath::Finite(h->GetXaxis()->GetBinUpEdge(b))) return false;
      return h->GetNbinsX() > 0 && h->GetXaxis()->GetXmin() >= 0;
    }
    if (g && g->GetN() >= 2) {
      double previous = -1;
      for (int j = 0; j < g->GetN(); ++j) {
        double x, y; g->GetPoint(j, x, y);
        if (!TMath::Finite(x) || !TMath::Finite(y) || x < 0 ||
            x <= previous || y < 0) return false;
        previous = x;
      }
      return true;
    }
    return false;
  }
  static bool evaluate(TObject *o, double x, double &value) {
    TH1 *h = dynamic_cast<TH1*>(o);
    TGraph *g = dynamic_cast<TGraph*>(o);
    if (h) {
      int b = h->FindFixBin(x);
      if (b < 1 || b > h->GetNbinsX()) return false;
      value = h->GetBinContent(b);
    } else if (g) {
      if (x < g->GetX()[0] || x > g->GetX()[g->GetN()-1]) return false;
      value = g->Eval(x);
    } else return false;
    return TMath::Finite(value) && value >= 0;
  }
  // Integrate the supplied support: piecewise constant TH1 or linear TGraph.
  static double integral(TObject *o, double upper) {
    double sum = 0;
    TH1 *h = dynamic_cast<TH1*>(o);
    TGraph *g = dynamic_cast<TGraph*>(o);
    if (h) for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double lo = h->GetXaxis()->GetBinLowEdge(b);
      double hi = TMath::Min(upper, h->GetXaxis()->GetBinUpEdge(b));
      if (hi > lo) sum += (hi-lo)*h->GetBinContent(b);
    }
    if (g) for (int j = 1; j < g->GetN(); ++j) {
      double lo = g->GetX()[j-1], hi = TMath::Min(upper, g->GetX()[j]);
      if (hi <= lo) continue;
      double ylo = g->GetY()[j-1];
      double yhi = ylo + (g->GetY()[j]-ylo)*(hi-lo)/(g->GetX()[j]-lo);
      sum += (hi-lo)*(ylo+yhi)/2;
    }
    return sum;
  }
};

// Load the three-column table written by reweight_egamma.C into the exact
// logarithmic E_gamma bins used by this analysis.  Keeping a TH1 lookup avoids
// interpolating a quantity that was defined bin by bin.
bool loadNoonTextWeights(const char *path, TH1D &weightHistogram)
{
  ifstream input(path);
  if (!input.is_open()) {
    cerr << "Cannot open n00n weight file: " << path << endl;
    return false;
  }

  string header;
  getline(input, header);
  vector<bool> seen(weightHistogram.GetNbinsX() + 1, false);
  int bin = 0, loadedBins = 0;
  double energyMeV = 0, weight = 0;
  while (input >> bin >> energyMeV >> weight) {
    if (bin < 1 || bin > weightHistogram.GetNbinsX() || seen[bin]) {
      cerr << "Invalid or duplicate n00n weight bin " << bin << endl;
      return false;
    }
    const double expectedEnergy = weightHistogram.GetBinCenter(bin);
    const double tolerance = 1.e-5*TMath::Max(1.0, TMath::Abs(expectedEnergy));
    if (!TMath::Finite(energyMeV) ||
        TMath::Abs(energyMeV-expectedEnergy) > tolerance) {
      cerr << "n00n weight energy mismatch in bin " << bin
           << ": file=" << energyMeV
           << " MeV, expected=" << expectedEnergy << " MeV" << endl;
      return false;
    }
    if (!TMath::Finite(weight) || weight < 0) {
      cerr << "Invalid n00n weight in bin " << bin << endl;
      return false;
    }
    weightHistogram.SetBinContent(bin, weight);
    weightHistogram.SetBinError(bin, 0);
    seen[bin] = true;
    ++loadedBins;
  }

  if (!input.eof() || loadedBins != weightHistogram.GetNbinsX()) {
    cerr << "Expected " << weightHistogram.GetNbinsX()
         << " n00n weight bins but loaded " << loadedBins << endl;
    return false;
  }
  cout << "Loaded " << loadedBins << " exact-bin n00n weights from "
       << path << endl;
  return true;
}

// Check the serialized output, including under/overflow and every tree leaf.
bool finiteOutput(TDirectory *directory)
{
  bool ok = true;
  TIter next(directory->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TObject *object = key->ReadObj();
    TH1 *h = dynamic_cast<TH1*>(object);
    TTree *tree = dynamic_cast<TTree*>(object);
    TDirectory *sub = dynamic_cast<TDirectory*>(object);
    bool objectOK = true;
    if (h) for (int b = 0; b < (h->GetNbinsX()+2)*(h->GetDimension()>1 ? h->GetNbinsY()+2 : 1)*(h->GetDimension()>2 ? h->GetNbinsZ()+2 : 1); ++b)
      if (!TMath::Finite(h->GetBinContent(b)) || !TMath::Finite(h->GetBinError(b))) objectOK = false;
    if (h) {
      double stats[13] = {0}; h->GetStats(stats);
      for (int j = 0; j < 13; ++j) if (!TMath::Finite(stats[j])) objectOK = false;
      TProfile *profile = dynamic_cast<TProfile*>(h);
      if (profile) for (int b = 0; b <= h->GetNbinsX()+1; ++b)
        if (!TMath::Finite(profile->GetBinEntries(b))) objectOK = false;
    }
    TGraph *graph = dynamic_cast<TGraph*>(object);
    if (graph) for (int j = 0; j < graph->GetN(); ++j)
      if (!TMath::Finite(graph->GetX()[j]) || !TMath::Finite(graph->GetY()[j])) objectOK = false;
    if (tree) for (Long64_t e = 0; e < tree->GetEntries(); ++e) {
      if (tree->GetEntry(e) <= 0) { objectOK = false; break; }
      TIter leaves(tree->GetListOfLeaves());
      TLeaf *leaf;
      while ((leaf = (TLeaf*)leaves()))
        for (int v = 0; v < leaf->GetNdata(); ++v)
          if (!TMath::Finite(leaf->GetValue(v))) objectOK = false;
    }
    if (sub) objectOK = finiteOutput(sub);
    if (!objectOK) cerr << "Nonfinite output: " << directory->GetPath() << "/" << object->GetName() << endl;
    ok = objectOK && ok;
    delete object;
  }
  return ok;
}

// define histogramshPartEta
TH1F *hEventCounter;
TH1F *hPartEta[10];

int main(int argc, char **argv)
{
  if(argc != 4 && argc != 5 && argc != 8) {
    cerr << "Usage: " << argv[0]
         << " <simulation type> <file list> <output file>"
         << " [reweighting_factors.txt]"
         << " OR [noon.root weightObject fluxDensityObject coverageThresholdGeV]"
         << endl;
    return 1;
  }
  TString simType = argv[1];
  TString fileList = argv[2];
  TString outFileName = argv[3];
  
  const bool useTextWeights = argc == 5;
  const TString textWeightPath = useTextWeights ? argv[4] : "";
  const double tpcEfficiency = 0.95;
  const double tofEfficiency = 0.65;
  const UInt_t tpcRandomSeed = 12345;
  const UInt_t tofRandomSeed = 67890;
  // This unweighted STAR sample stores an unset scalar weight as zero.
  // Set false for samples in which zero is an intentional generator weight.
  const bool zeroGeneratorWeightIsUnset = true;
  TRandom3 tpcRandom(tpcRandomSeed);
  TRandom3 tofRandom(tofRandomSeed);
  NoonFluxInput flux;
  double coverageThresholdGeV = 0;
  double fluxFractionBelowCoverage = 0;
  if (argc == 8) {
    char *end = 0;
    coverageThresholdGeV = strtod(argv[7], &end);
    flux.file = TFile::Open(argv[4], "READ");
    if (!flux.file || flux.file->IsZombie() || !end || *end || end == argv[7] ||
        !TMath::Finite(coverageThresholdGeV) || coverageThresholdGeV < 0) {
      cerr << "Invalid nOOn input file or coverage threshold" << endl; return 1;
    }
    flux.weights = flux.file->Get(argv[5]);
    flux.density = flux.file->Get(argv[6]);
    if (!NoonFluxInput::validate(flux.weights) || !NoonFluxInput::validate(flux.density)) {
      cerr << "nOOn objects must be finite nonnegative 1D TH1 or strictly increasing TGraph (x in MeV)" << endl;
      return 1;
    }
    double allFlux = NoonFluxInput::integral(flux.density, TMath::Infinity());
    if (!TMath::Finite(allFlux) || allFlux <= 0) {
      cerr << "Invalid or zero integrated physical nOOn flux" << endl; return 1;
    }
    fluxFractionBelowCoverage = NoonFluxInput::integral(flux.density, coverageThresholdGeV*1000)/allFlux;
  } else if (!useTextWeights) {
    cout << "No nOOn flux input: output is unweighted by UPC flux "
         << "(noonFluxWeight=1)." << endl;
  }

  //----------------------------------
  // Open files and add to the chain
  //----------------------------------
  printf("\n----------------------------------\n");
  printf("Open files and add to the chain\n");
  printf("----------------------------------\n\n");
  TChain *chain = NULL;
  chain = new TChain("genevents");

  int ifile=0;
  char filename[512];
  ifstream *inputStream = new ifstream(fileList.Data());
  if (!inputStream->is_open()) 
    {
      printf("Cannot open file list: %s\n",fileList.Data());
      return 0;
    }
  for(;inputStream->good();)
    {
      inputStream->getline(filename,512);
      if(inputStream->good()) 
	{
	  if (strlen(filename) == 0) continue;
	  TFile *ftmp = TFile::Open(filename, "READ");
	  if(!ftmp || ftmp->IsZombie() || !(ftmp->GetNkeys())) 
	    {
	      cout<<"something wrong"<<endl;
	    } 
	  else
	    {
	      cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
	      chain->Add(filename);
	      ifile++;
	    }
	  if (ftmp) {
	    ftmp->Close();
	    delete ftmp;
	  }
	}
    }
  delete inputStream;

  //------------------------------------------------
  // Define histograms
  //------------------------------------------------
  printf("\n----------------------------------\n");
  printf("Define histograms\n");
  printf("----------------------------------\n\n");

  // Open before creating the event tree so baskets stream to disk.
  TFile *fout = new TFile(outFileName.Data(), "RECREATE");
  if (fout->IsZombie()) { cerr << "Cannot create output" << endl; return 1; }
  fout->cd();
  TH1::SetDefaultSumw2();
  hEventCounter = new TH1F("hEventCounter", "Event statistics;sub-process", 110, 0, 110);

  const int nBins = 300;
  double logBins[nBins + 1];
  const double log_kMin = log10(10.0), log_kMax = log10(5e6);

  for (int i = 0; i <= nBins; ++i) {
         logBins[i]  = pow(10, log_kMin + (log_kMax - log_kMin) * i / nBins);
  }

  TH1D hNoonTextWeight(
      "hNoonTextWeight",
      "Applied n00n exact-bin weight;E_{#gamma} [MeV];Weight",
      nBins, logBins);
  if (useTextWeights &&
      !loadNoonTextWeights(textWeightPath.Data(), hNoonTextWeight)) return 1;
  TObject *activeWeightObject =
      useTextWeights ? static_cast<TObject*>(&hNoonTextWeight) : flux.weights;

  TH1D fullq2hist("fullq2hist","Q^2 [GeV^2]", 100, 0.,0.1);
  TH1D egamma("egamma","E_{#gamma} [MeV]", nBins, logBins);
  TH1D q2hist("q2hist","Q^2 [GeV^2]", 100, 0., 0.01);
  TH1D ypionplus("ypionplus","",120,-6.,6.);
  TH1D ykaonplus("ykaonplus","",120,-6.,6.);
  TH1D yprotonplus("yprotonplus","",120,-6.,6.);
  TH1D ypionmins("ypionmins","",120,-6.,6.);
  TH1D ykaonmins("ykaonmins","",120,-6.,6.);
  TH1D yprotonmins("yprotonmins","",120,-6.,6.);
  TH1D yneutron("yneutron","",120,-6.,6.);


  TH1D ptpionplus("ptpionplus","",200,0.,2.);
  TH1D ptkaonplus("ptkaonplus","",200,0.,2.);
  TH1D ptprotonplus("ptprotonplus","",200,0.,2.);
  TH1D ptpionmins("ptpionmins","",200,0.,2.);
  TH1D ptkaonmins("ptkaonmins","",200,0.,2.);
  TH1D ptprotonmins("ptprotonmins","",200,0.,2.);
  TH1D ptneutron("ptneutron","",200,0.,2.);

  TH1D *h_dNdEta = new TH1D("h_dNdEta", "Raw charged final-state particle yield;#eta;Weighted yield", 201, -5.025, 5.025);
  TH1D *h_dNdEta2 = new TH1D("h_dNdEta2", "TPC-reconstructed charged tracks, |#eta|<1, 0.2<p_{T}<2.0;#eta;Weighted yield", 201, -5.025, 5.025);
  TH1D hAcceptedEventWeight("hAcceptedEventWeight", "Accepted event weight sum", 1, 0, 1);
  TH1D hAcceptedEvents("hAcceptedEvents", "Accepted event count", 1, 0, 1);
  TH1D hPhotonUnresolved("hPhotonUnresolved", "Missing or ambiguous beam photon", 1, 0, 1);
  const int speciesPdg[7] = {211, -211, 321, -321, 2212, -2212, 2112};
  const char *speciesName[7] = {"pionplus", "pionminus", "kaonplus", "kaonminus", "proton", "antiproton", "neutron"};
  TH1D *hEtaSpecies[7];
  TH1D *hShiftedY[6];
  const double Ybeam = acosh(27.0 / 0.938272);
  for (int species = 0; species < 7; ++species) {
    hEtaSpecies[species] = new TH1D(Form("hEta_%s", speciesName[species]),
        "Inclusive species, p_{T}>0.2;#eta;Weighted yield", 160, -8, 8);
    if (species < 6) hShiftedY[species] = new TH1D(Form("hShiftedY_%s", speciesName[species]),
        "0.4<p_{T}<1.2;y-Y_{beam};Weighted yield", 160, -12, 4);
  }

  TH1D *h_pt_all = new TH1D("h_pt_all", "TPC-reconstructed charged tracks;p_{T} [GeV/c];Weighted yield", 200, 0., 5.);

  TH1D *h_npart_raw = new TH1D("h_npart_raw", "Raw Final-State Multiplicity;Multiplicity;Events", 200, 0, 200);
  TH1D *h_npart_tof_matched = new TH1D("h_npart_tof_matched", "TOF-Matched Multiplicity;Multiplicity;Events", 200, -0.5, 199.5);

  TH2D *h_deta_dphi = new TH2D("h_deta_dphi", "#Delta#eta-#Delta#phi Correlations;#Delta#eta;#Delta#phi", 60, 0.0, 3.0, 24, -TMath::Pi()/2, 3*TMath::Pi()/2);
  const int nCentBins = 3;
  TH1D* h_dNdEta_cent[nCentBins];
  TH1D* h_dNdEta2_cent[nCentBins];
  TH2D* h_deta_dphi_cent[nCentBins];  
  const char* centBinTitles[] = {"(N_{ToF} <= 5)", "(5 < N_{ToF} <= 15)", "(N_{ToF} > 15)"};
  
  for (int i = 0; i < nCentBins; ++i) {
    h_dNdEta_cent[i] = new TH1D(
        Form("h_dNdEta_cent%d", i),
        Form("Raw charged final-state particle yield %s;#eta;Weighted yield", centBinTitles[i]),
        201, -5.025, 5.025);
    h_dNdEta2_cent[i] = new TH1D(
        Form("h_dNdEta2_cent%d", i),
        Form("TPC-reconstructed charged tracks, |#eta|<1, 0.2<p_{T}<2.0 %s;#eta;Weighted yield", centBinTitles[i]),
        201, -5.025, 5.025);
    h_deta_dphi_cent[i] = new TH2D(
        Form("h_deta_dphi_cent%d", i), 
	Form("#Delta#eta-#Delta#phi %s;#Delta#eta;#Delta#phi", centBinTitles[i]), 60, 0.0, 3.0, 24, -TMath::Pi()/2, 3*TMath::Pi()/2);
  }  

  const char* partName[10] = {"charge", "baryon", "pion_plus", "kaon_plus", "proton_plus", "pion_minus", "kaon_minus", "proton_minus", "proton_plus_prim", "proton_minus_prim"};
  for(int i=0; i<10; i++)
  {
	  hPartEta[i] = new TH1F(Form("hPartEta_%s",partName[i]), "Particle rapidity distribution; y", 200, -10, 10);
	  hPartEta[i]->Sumw2();
  }

  // Preserve the legacy graph output name for plotting/QA.  Text weights are
  // still evaluated through hNoonTextWeight, not by graph interpolation.
  TGraph *gWeights = new TGraph();
  gWeights->SetName("gWeights");
  gWeights->SetTitle("Applied UPC weight;E_{#gamma} [MeV];Weight");
  if (!activeWeightObject) {
    gWeights->SetPoint(0, 0, 1); gWeights->SetPoint(1, 5000, 1);
  } else if (TGraph *g = dynamic_cast<TGraph*>(activeWeightObject)) {
    for (int j = 0; j < g->GetN(); ++j) gWeights->SetPoint(j, g->GetX()[j], g->GetY()[j]);
  } else {
    TH1 *h = (TH1*)activeWeightObject;
    for (int j = 1; j <= h->GetNbinsX(); ++j) gWeights->SetPoint(j-1, h->GetBinCenter(j), h->GetBinContent(j));
  }

  TH2D h2_Egamma_Ntof_generator("h2_Egamma_Ntof_generator", ";E_{#gamma} [MeV];N_{trk}^{TOF}", nBins, logBins, 200, -0.5, 199.5);
  TH2D h2_Egamma_Ntof_upc("h2_Egamma_Ntof_upc", ";E_{#gamma} [MeV];N_{trk}^{TOF}", nBins, logBins, 200, -0.5, 199.5);
  TH1D h_Egamma_generator_all("h_Egamma_generator_all", ";E_{#gamma} [MeV];Generator yield", nBins, logBins);
  TH1D h_Egamma_upc_all("h_Egamma_upc_all", ";E_{#gamma} [MeV];UPC yield", nBins, logBins);
  TH1D *generatorClass[3], *upcClass[3], *probability[3];
  for (int c = 0; c < 3; ++c) {
    generatorClass[c] = new TH1D(Form("h_Egamma_generator_class%d", c), ";E_{#gamma} [MeV];Generator yield", nBins, logBins);
    upcClass[c] = new TH1D(Form("h_Egamma_upc_class%d", c), ";E_{#gamma} [MeV];UPC yield", nBins, logBins);
    probability[c] = new TH1D(Form("h_class%d_probability_vs_Egamma", c), ";E_{#gamma} [MeV];P(class|E_{#gamma})", nBins, logBins);
  }
  TProfile p_meanNtof_vs_Egamma("p_meanNtof_vs_Egamma", ";E_{#gamma} [MeV];Mean N_{trk}^{TOF} (generator weight)", nBins, logBins);
  TH1D h_Ntof_generator("h_Ntof_generator", ";N_{trk}^{TOF};Generator yield", 200, -0.5, 199.5);
  TH1D h_Ntof_upc("h_Ntof_upc", ";N_{trk}^{TOF};n00n-weighted yield", 200, -0.5, 199.5);

  TH1D h_Nraw_generator("h_Nraw_generator", "Raw charged multiplicity;N_{ch}^{raw};Generator yield", 200, -0.5, 199.5);
  TH1D h_Nraw_upc("h_Nraw_upc", "Raw charged multiplicity;N_{ch}^{raw};n00n-weighted yield", 200, -0.5, 199.5);
  TH1D h_Naccept_generator("h_Naccept_generator", "Acceptance-selected multiplicity;N_{ch}^{accept};Generator yield", 200, -0.5, 199.5);
  TH1D h_Naccept_upc("h_Naccept_upc", "Acceptance-selected multiplicity;N_{ch}^{accept};n00n-weighted yield", 200, -0.5, 199.5);
  TH1D h_Ntpc_generator("h_Ntpc_generator", "TPC-reconstructed multiplicity;N_{trk}^{TPC};Generator yield", 200, -0.5, 199.5);
  TH1D h_Ntpc_upc("h_Ntpc_upc", "TPC-reconstructed multiplicity;N_{trk}^{TPC};n00n-weighted yield", 200, -0.5, 199.5);

  TH2D h2_Egamma_Nraw_generator("h2_Egamma_Nraw_generator", ";E_{#gamma} [MeV];N_{ch}^{raw}", nBins, logBins, 200, -0.5, 199.5);
  TH2D h2_Egamma_Nraw_upc("h2_Egamma_Nraw_upc", ";E_{#gamma} [MeV];N_{ch}^{raw}", nBins, logBins, 200, -0.5, 199.5);
  TH2D h2_Egamma_Naccept_generator("h2_Egamma_Naccept_generator", ";E_{#gamma} [MeV];N_{ch}^{accept}", nBins, logBins, 200, -0.5, 199.5);
  TH2D h2_Egamma_Naccept_upc("h2_Egamma_Naccept_upc", ";E_{#gamma} [MeV];N_{ch}^{accept}", nBins, logBins, 200, -0.5, 199.5);
  TH2D h2_Egamma_Ntpc_generator("h2_Egamma_Ntpc_generator", ";E_{#gamma} [MeV];N_{trk}^{TPC}", nBins, logBins, 200, -0.5, 199.5);
  TH2D h2_Egamma_Ntpc_upc("h2_Egamma_Ntpc_upc", ";E_{#gamma} [MeV];N_{trk}^{TPC}", nBins, logBins, 200, -0.5, 199.5);

  const char *responseStageName[4] = {"raw", "accept", "tpc", "tof"};
  const char *responseStageTitle[4] = {
      "at least one raw charged particle",
      "at least one acceptance-selected charged particle",
      "at least one TPC-reconstructed track",
      "at least one TOF-matched track"};
  TH1D *h_Egamma_stage_generator[4];
  TH1D *h_Egamma_stage_upc[4];
  for (int stage = 0; stage < 4; ++stage) {
    h_Egamma_stage_generator[stage] = new TH1D(
        Form("h_Egamma_%s_generator", responseStageName[stage]),
        Form("E_{#gamma}, %s;E_{#gamma} [MeV];Generator yield", responseStageTitle[stage]),
        nBins, logBins);
    h_Egamma_stage_upc[stage] = new TH1D(
        Form("h_Egamma_%s_upc", responseStageName[stage]),
        Form("E_{#gamma}, %s;E_{#gamma} [MeV];n00n-weighted yield", responseStageTitle[stage]),
        nBins, logBins);
  }
  TH1D h_WgammaP_all("h_WgammaP_all", ";W_{#gamma p} [GeV];UPC yield", 400, 0, 200);
  TH1D h_WgammaP_class0("h_WgammaP_class0", ";W_{#gamma p} [GeV];UPC yield", 400, 0, 200);
  TH1D h_WgammaP_class1("h_WgammaP_class1", ";W_{#gamma p} [GeV];UPC yield", 400, 0, 200);
  TH1D hProtonUnresolved("hProtonUnresolved", "Missing, ambiguous or invalid incoming proton", 1, 0, 1);
  Long64_t eventIndex = 0;
  double egammaGeV = 0, egammaMeV = 0, q2GeV2 = 0, wGammaPGeV = 0;
  int photonDirection = 0, nFinalState = 0, nRawCharged = 0;
  int nTracksAccepted = 0, nTracksTpc = 0, nTofMatched = 0, multClass = 0;
  double generatorWeight = 1, noonFluxWeight = 1, totalWeight = 1;
  TTree noonEvents("noonEvents", "nOOn interface: GenerateEvent(egammaGeV), separately for multClass 0 and 1; class 2 retained for closure");
  noonEvents.Branch("eventIndex", &eventIndex, "eventIndex/L");
  noonEvents.Branch("egammaGeV", &egammaGeV, "egammaGeV/D");
  noonEvents.Branch("egammaMeV", &egammaMeV, "egammaMeV/D");
  noonEvents.Branch("q2GeV2", &q2GeV2, "q2GeV2/D");
  noonEvents.Branch("wGammaPGeV", &wGammaPGeV, "wGammaPGeV/D");
  noonEvents.Branch("photonDirection", &photonDirection, "photonDirection/I");
  noonEvents.Branch("nFinalState", &nFinalState, "nFinalState/I");
  noonEvents.Branch("nRawCharged", &nRawCharged, "nRawCharged/I");
  noonEvents.Branch("nTracksAccepted", &nTracksAccepted, "nTracksAccepted/I");
  noonEvents.Branch("nTracksTpc", &nTracksTpc, "nTracksTpc/I");
  noonEvents.Branch("nTofMatched", &nTofMatched, "nTofMatched/I");
  noonEvents.Branch("multClass", &multClass, "multClass/I");
  noonEvents.Branch("generatorWeight", &generatorWeight, "generatorWeight/D");
  noonEvents.Branch("noonFluxWeight", &noonFluxWeight, "noonFluxWeight/D");
  noonEvents.Branch("totalWeight", &totalWeight, "totalWeight/D");
  Long64_t classCount[3] = {0, 0, 0}, invalidGeneratorWeights = 0;
  double classWeight[3] = {0, 0, 0};
  bool unitGeneratorWeights = true;

  //------------------------------------------------
  // Loop the chain
  //------------------------------------------------
  Long64_t nevents = chain->GetEntries();
  printf("\n----------------------------------\n");
  printf("Loop the chain: %lld events\n",(long long)nevents);
  printf("----------------------------------\n\n");

  genevents *t = new genevents(chain);

  for(Long64_t i=0;i<nevents;i++)
  {
	  //if(i==2000) break;
	  if(i%1000==0) cout << "begin " << i << "th entry...." << endl;
	  if (t->LoadTree(i) < 0) { cerr << "Cannot load event " << i << endl; return 1; }
          // Read the scalar generator weight from the current underlying tree,
          // avoiding stale values if a file in the chain lacks this branch.
          t->weight = 1.0;
          if (t->GetEntry(i) <= 0) { cerr << "Cannot read event " << i << endl; return 1; }
          generatorWeight = 1.0;
          TLeaf *weightLeaf = chain->GetTree()->GetLeaf("weight");
          if (weightLeaf && weightLeaf->GetLen() == 1) {
            double candidate = weightLeaf->GetValue();
            if (TMath::Finite(candidate) && !(zeroGeneratorWeightIsUnset && candidate == 0)) generatorWeight = candidate;
            else ++invalidGeneratorWeights;
          }
          unitGeneratorWeights = unitGeneratorWeights && generatorWeight == 1.0;

	  /*
	  // remove elastic process
	  if(simType.Contains("py6"))
	  {
	  if(t->subprocess==91) continue;
	  else hEventCounter->Fill(t->subprocess);
	  }
	  if(simType.Contains("py8"))
	  {
	  if(t->process==102) continue;
	  else hEventCounter->Fill(t->process);
	  }
	  */

	  int n_forward_tracks = 0;
	  int n_backward_tracks = 0;
	  std::vector<TLorentzVector> forward_tracks;
	  std::vector<TLorentzVector> backward_tracks;

	  int npart = t->mNumParticles;
	  //printf("\n Found %d particles \n",npart);
	  int net_charge = 0;
	  int net_baryon = 0;


          TLorentzVector pIncoming;
          const bool hasIncomingProton = incomingProton(t, pIncoming);
          if (!hasIncomingProton) hProtonUnresolved.Fill(0.5);
          TLorentzVector q;
          if (!exchangedPhoton(t, q) || !TMath::Finite(q.E()) ||
              !TMath::Finite(q.Px()) || !TMath::Finite(q.Py()) || !TMath::Finite(q.Pz()) || !TMath::Finite(q.M2())) {
            hPhotonUnresolved.Fill(0.5);
            continue;
          }
          const double E_gamma = q.E(); // beam/laboratory frame
          const double Q2_gamma = -q.M2(); // Q2fac is NOT photon virtuality
          fullq2hist.Fill(Q2_gamma, generatorWeight);
          
	  if (!(E_gamma > 0.0 && E_gamma < 5.0)) continue;
          if (!(Q2_gamma >= 0.0 && Q2_gamma < 0.01)) continue;
          
	  if (!hasIncomingProton) continue;
          eventIndex = i;
          egammaGeV = E_gamma;
          egammaMeV = E_gamma*1000;
          q2GeV2 = Q2_gamma;
          wGammaPGeV = (q + pIncoming).M();
          if (!TMath::Finite(wGammaPGeV)) {
            cerr << "Invalid W_gamma_p at event " << i << endl; return 1;
          }
          photonDirection = (q.Pz() > 0) - (q.Pz() < 0);
          noonFluxWeight = 1.0;
          if (activeWeightObject &&
              !NoonFluxInput::evaluate(activeWeightObject, egammaMeV, noonFluxWeight)) {
            cerr << "nOOn weight lacks valid coverage at " << egammaMeV << " MeV; aborting" << endl; return 1;
          }
          totalWeight = generatorWeight*noonFluxWeight;
          if (!TMath::Finite(totalWeight)) { cerr << "Nonfinite event weight" << endl; return 1; }
          // UPC physics outputs use the generator weight multiplied by the
          // exact-bin n00n flux-shape weight. Generator-reference outputs keep
          // generatorWeight explicitly below.
          const double event_weight = totalWeight;
          hAcceptedEventWeight.Fill(0.5, event_weight);
          hAcceptedEvents.Fill(0.5);
          egamma.Fill(E_gamma*1e3, event_weight);
          q2hist.Fill(Q2_gamma, event_weight);

          int final_state_npart = 0;
	  for (int j = 0; j < npart; j++) {
    	    if (t->mParticles_mStatus[j] == 1) {
              final_state_npart++;
    	    }
	  }

	  h_npart_raw->Fill(final_state_npart, event_weight);

	  // Detector response chain:
	  // raw charged -> kinematic acceptance -> TPC reconstruction -> TOF match.
	  const double track_eta_max = 1.0;
	  const double track_pt_min = 0.2;
	  const double track_pt_max = 2.0;
	  std::vector<TLorentzVector> accepted_truth_tracks;
	  std::vector<TLorentzVector> tpc_tracks;
	  nRawCharged = 0;

	  for (int j = 0; j < npart; ++j) {
	    if (t->mParticles_mStatus[j] != 1) continue;
	    if (!isChargedParticle(t->mParticles_mId[j])) continue;
	    ++nRawCharged;

	    TLorentzVector track;
	    track.SetPxPyPzE(t->mParticles_mPx[j], t->mParticles_mPy[j],
	                    t->mParticles_mPz[j], t->mParticles_mEnergy[j]);

	    if (!TMath::Finite(track.Eta()) || !TMath::Finite(track.Pt())) continue;
	    if (fabs(track.Eta()) >= track_eta_max) continue;
	    if (track.Pt() <= track_pt_min || track.Pt() >= track_pt_max) continue;
	    accepted_truth_tracks.push_back(track);

	    // Independent constant TPC response. Replace this Bernoulli probability
	    // with a validated efficiency map epsilon_TPC(pT,eta) when available.
	    if (tpcRandom.Rndm() < tpcEfficiency) tpc_tracks.push_back(track);
	  }

	  nTracksAccepted = static_cast<int>(accepted_truth_tracks.size());
	  nTracksTpc = static_cast<int>(tpc_tracks.size());
          nTofMatched = 0;
          // TOF matching is conditional on successful TPC reconstruction.
          for (size_t k = 0; k < tpc_tracks.size(); ++k)
            if (tofRandom.Rndm() < tofEfficiency) ++nTofMatched;
          int npart_tof_matched = nTofMatched;
	  h_npart_tof_matched->Fill(npart_tof_matched, event_weight);

	  h_Nraw_generator.Fill(nRawCharged, generatorWeight);
	  h_Nraw_upc.Fill(nRawCharged, totalWeight);
	  h_Naccept_generator.Fill(nTracksAccepted, generatorWeight);
	  h_Naccept_upc.Fill(nTracksAccepted, totalWeight);
	  h_Ntpc_generator.Fill(nTracksTpc, generatorWeight);
	  h_Ntpc_upc.Fill(nTracksTpc, totalWeight);
	  h_Ntof_generator.Fill(nTofMatched, generatorWeight);
	  h_Ntof_upc.Fill(nTofMatched, totalWeight);

	  h2_Egamma_Nraw_generator.Fill(egammaMeV, nRawCharged, generatorWeight);
	  h2_Egamma_Nraw_upc.Fill(egammaMeV, nRawCharged, totalWeight);
	  h2_Egamma_Naccept_generator.Fill(egammaMeV, nTracksAccepted, generatorWeight);
	  h2_Egamma_Naccept_upc.Fill(egammaMeV, nTracksAccepted, totalWeight);
	  h2_Egamma_Ntpc_generator.Fill(egammaMeV, nTracksTpc, generatorWeight);
	  h2_Egamma_Ntpc_upc.Fill(egammaMeV, nTracksTpc, totalWeight);

	  const int responseCounts[4] = {
	      nRawCharged, nTracksAccepted, nTracksTpc, nTofMatched};
	  for (int stage = 0; stage < 4; ++stage) {
	    if (responseCounts[stage] <= 0) continue;
	    h_Egamma_stage_generator[stage]->Fill(egammaMeV, generatorWeight);
	    h_Egamma_stage_upc[stage]->Fill(egammaMeV, totalWeight);
	  }

	  int icent = -1;
          if (npart_tof_matched <= 5) {
            icent = 0;
	  } else if (npart_tof_matched > 5 && npart_tof_matched <= 15) {
            icent = 1;
	  } else {
	    icent = 2;
	  }


          multClass = icent;
          nFinalState = final_state_npart;
          ++classCount[multClass];
          classWeight[multClass] += totalWeight;
          h2_Egamma_Ntof_generator.Fill(egammaMeV, nTofMatched, generatorWeight);
          h2_Egamma_Ntof_upc.Fill(egammaMeV, nTofMatched, totalWeight);
          h_Egamma_generator_all.Fill(egammaMeV, generatorWeight);
          generatorClass[multClass]->Fill(egammaMeV, generatorWeight);
          h_Egamma_upc_all.Fill(egammaMeV, totalWeight);
          upcClass[multClass]->Fill(egammaMeV, totalWeight);
          p_meanNtof_vs_Egamma.Fill(egammaMeV, nTofMatched, generatorWeight);
          h_WgammaP_all.Fill(wGammaPGeV, totalWeight);
          if (multClass == 0) h_WgammaP_class0.Fill(wGammaPGeV, totalWeight);
          if (multClass == 1) h_WgammaP_class1.Fill(wGammaPGeV, totalWeight);
          if (noonEvents.Fill() < 0) { cerr << "Cannot write event tree" << endl; return 1; }

          // Analysis tracks are TPC-reconstructed tracks. TOF matching defines
          // the event-activity class but does not further thin pair particles.
          for (const auto& track : tpc_tracks) {
            h_dNdEta2->Fill(track.Eta(), event_weight);
            h_dNdEta2_cent[icent]->Fill(track.Eta(), event_weight);
            h_pt_all->Fill(track.Pt(), event_weight);
          }


	  /*for (int j = 0; j < npart_tof_matched; j++) {
            if (t->mParticles_mStatus[j] != 1) continue;

            TLorentzVector mom;
            mom.SetPxPyPzE(t->mParticles_mPx[j], t->mParticles_mPy[j], t->mParticles_mPz[j], t->mParticles_mEnergy[j]);
            double eta = mom.Eta();
	    //cout<<"j_npart="<<j<<" eta="<<eta<<endl;
            if (eta > 2.1 && eta < 5.0) {
                n_forward_tracks++;
                forward_tracks.push_back(mom);
            }
            else if (eta > -5.0 && eta < -2.1) {
                n_backward_tracks++;
                backward_tracks.push_back(mom);
            }
          }
	  //cout << "nforward=" << n_forward_tracks << " " << "nbackward=" << n_backward_tracks << endl;

          std::vector<TLorentzVector> selected_tracks;
          if (n_forward_tracks > 1 && n_backward_tracks <= 1) {
              selected_tracks = forward_tracks;
          }
          else if (n_backward_tracks > 1 && n_forward_tracks <= 1) {
              selected_tracks = backward_tracks;
          }*/

	  // Fill inclusive spectra from every final-state particle.
	  for(int j=0; j<npart; j++)
	  {
		  if(t->mParticles_mStatus[j]!=1) continue;
		  int pdg = t->mParticles_mId[j];
		  TLorentzVector mom;
		  mom.SetPxPyPzE(t->mParticles_mPx[j], t->mParticles_mPy[j], t->mParticles_mPz[j], t->mParticles_mEnergy[j]);
		  double rap = mom.Rapidity(); //change of convention, nothing changes
		  double pt = mom.Pt();

          // All charged final-state particles, without track acceptance cuts.
          if (isChargedParticle(pdg)) {
            h_dNdEta->Fill(mom.Eta(), event_weight);
            h_dNdEta_cent[icent]->Fill(mom.Eta(), event_weight);
          }
          const double dy = mom.Rapidity() - Ybeam;
          for (int species = 0; species < 7; ++species) {
            if (pdg != speciesPdg[species]) continue;
            if (pt > 0.2) hEtaSpecies[species]->Fill(mom.Eta(), event_weight);
            if (species < 6 && pt > 0.2 && pt < 2.0)// pt range changed from 0.4 to 0.2, and from 1.2 to 2.0
              hShiftedY[species]->Fill(dy, event_weight);
          }
		  //cout<<" id = 2212 "<<pdg<<" y= "<<rap<<" pt= "<<pt<<endl;

		  //if(pt>0.4 && pt<1.2){
		  if(pdg==211)ypionplus.Fill(rap, event_weight);
		  if(pdg==321)ykaonplus.Fill(rap, event_weight);
		  if(pdg==2212)yprotonplus.Fill(rap, event_weight);
		  if(pdg==-211)ypionmins.Fill(rap, event_weight);
		  if(pdg==-321)ykaonmins.Fill(rap, event_weight);
		  if(pdg==-2212)yprotonmins.Fill(rap, event_weight);
		  if(pdg==2112)yneutron.Fill(rap, event_weight);

		  if(pdg==211)ptpionplus.Fill(pt, event_weight);
		  if(pdg==321)ptkaonplus.Fill(pt, event_weight);
		  if(pdg==2212)ptprotonplus.Fill(pt, event_weight);

		  if(pdg==-211)ptpionmins.Fill(pt, event_weight);
		  if(pdg==-321)ptkaonmins.Fill(pt, event_weight);
		  if(pdg==-2212)ptprotonmins.Fill(pt, event_weight);

		  if(pdg==2112)ptneutron.Fill(pt, event_weight);


		  //}


		  //if(!isPrimaryParticle(pdg)) continue;
		  int charge = 0;      
		  int baryon = 0;
		  if(pdg==11)         { charge = -1; baryon = 0;  }
		  if(pdg==-11)        { charge = 1;  baryon = 0;  }
		  if(pdg==13)         { charge = -1; baryon = 0;  }
		  if(pdg==-13)        { charge = 1;  baryon = 0;  }
		  if(fabs(pdg)==12)   { charge = 0;  baryon = 0;  }
		  if(fabs(pdg)==14)   { charge = 0;  baryon = 0;  }
		  if(fabs(pdg)==16)   { charge = 0;  baryon = 0;  }
		  if(pdg==22)         { charge = 0;  baryon = 0;  }
		  if(pdg==111)        { charge = 0;  baryon = 0;  }
		  if(pdg==130)        { charge = 0;  baryon = 0;  }
		  if(pdg==211)        { charge = 1;  baryon = 0;  }
		  if(pdg==-211)       { charge = -1; baryon = 0;  }
		  if(pdg==310)        { charge = 0;  baryon = 0;  }
		  if(pdg==321)        { charge = 1;  baryon = 0;  }
		  if(pdg==-321)       { charge = -1; baryon = 0;  }
		  if(pdg==2112)       { charge = 0;  baryon = 1;  }
		  if(pdg==-2112)      { charge = 0;  baryon = -1; }
		  if(pdg==2212)       { charge = 1;  baryon = 1;  }
		  if(pdg==-2212)      { charge = -1; baryon = -1; }
		  // if(pdg==3112)       { charge = -1; baryon = 1;  }
		  // if(pdg==-3112)      { charge = 1; baryon = -1;  }
		  // if(pdg==3212)       { charge = 0; baryon = 1;   }
		  // if(pdg==-3212)      { charge = 0; baryon = -1;  }
		  // if(pdg==3222)       { charge = 1; baryon = 1;   }
		  // if(pdg==-3222)      { charge = -1; baryon = -1; }
		  // if(pdg==3122)       { charge = 0; baryon = 1;   }
		  // if(pdg==-3122)      { charge = 0; baryon = -1;  }
		  // if(pdg==3322)       { charge = 0; baryon = 1;   }
		  // if(pdg==-3322)      { charge = 0; baryon = -1;  }
		  // if(pdg==3312)       { charge = -1; baryon = 1;  }
		  // if(pdg==-3312)      { charge = 1; baryon = -1;  }
		  // if(pdg==3334)       { charge = -1; baryon = 1;  }
		  // if(pdg==-3334)      { charge = 1; baryon = -1;  }

		  //if(fabs(pdg==211)) cout << "Found a pion" << endl;

		  net_charge += charge;
		  net_baryon += baryon;

		  // if(fabs(pdg)==3122)
		  //   {
		  //     cout << "Found a Lambda: " << pdg << endl;
		  //     int daug_index_1 = t->mParticles_mDaughter[j][0];
		  //     int daug_index_2 = t->mParticles_mDaughter[j][1];
		  //     if(daug_index_1>=0)
		  //     //if(mom_index_1>=0)
		  // 	{
		  // 	  cout << "Daughter 1: " << t->mParticles_mStatus[daug_index_1] << "  " << t->mParticles_mId[daug_index_1] << endl;
		  // 	}
		  //     if(daug_index_2>=0)
		  // 	{
		  // 	  cout << "Daughter 2: " << t->mParticles_mStatus[daug_index_2] << "  " << t->mParticles_mId[daug_index_2] << endl;
		  // 	}
		  //   }

		  //if(baryon==-999)
		  //cout << j << ": pdg = " << pdg << ", charge = " << charge << ", baryon = " << baryon << endl;

		  //if(fabs(pdg)>1000)cout << j << ": pdg = " << pdg << ", mother = (" << t->mParticles_mMother[j][0] << ", " << t->mParticles_mMother[j][1] << "), daughter = (" << t->mParticles_mDaughter[j][0] << ", " << t->mParticles_mDaughter[j][1] << ")" << endl;

		  if(fabs(pdg)<100) continue;

		  hPartEta[0]->Fill(rap, charge);
		  hPartEta[1]->Fill(rap, baryon);

		  if(fabs(pdg)==211)
		  {
			  if(pdg==211)  hPartEta[2]->Fill(rap);
			  if(pdg==-211) hPartEta[5]->Fill(rap);
		  }

		  if(fabs(pdg)==321)
		  {
			  if(pdg==321)  hPartEta[3]->Fill(rap);
			  if(pdg==-321) hPartEta[6]->Fill(rap);
		  }

		  if(fabs(pdg)==2212)
		  {
			  if(pdg==2212)  hPartEta[4]->Fill(rap);
			  if(pdg==-2212) hPartEta[7]->Fill(rap);

			  bool is_weak = false;
			  int mom_index_1 = t->mParticles_mMother[j][0];
			  if(mom_index_1>=0)
			  {
				  int pdg_mom = t->mParticles_mId[mom_index_1];
				  if(fabs(pdg_mom)==3112 || fabs(pdg_mom)==3212 || fabs(pdg_mom)==3222 || fabs(pdg_mom)==3122 || fabs(pdg_mom)==3322 || fabs(pdg_mom)==3312 || fabs(pdg_mom)==3334)
					  is_weak = true;
			  }
			  if(is_weak==false)
			  {
				  if(pdg==2212)  hPartEta[8]->Fill(rap);
				  if(pdg==-2212) hPartEta[9]->Fill(rap);
			  }
		  }
		  //cout << j << ": " << t->mParticles_mStatus[j] << "  " << t->mParticles_mId[j] << "  " << t->mParticles_mMother[j][0] << "  " <<
		  // "  " << t->mParticles_mMass[j] << endl;
	  }

	  // Evaluate each physical pair once and fill both trigger/associate
	  // orientations to retain the symmetric correlation density.
	  for (size_t j = 0; j < tpc_tracks.size(); ++j) {
	    for (size_t k = j + 1; k < tpc_tracks.size(); ++k) {
	      double deta = tpc_tracks[k].Eta() - tpc_tracks[j].Eta();
	      double dphi = TVector2::Phi_mpi_pi(
	          tpc_tracks[k].Phi() - tpc_tracks[j].Phi());

	      if (dphi < -TMath::Pi()/2) dphi += 2*TMath::Pi();
	      h_deta_dphi->Fill(deta, dphi, event_weight);
	      h_deta_dphi_cent[icent]->Fill(deta, dphi, event_weight);

	      double reverse_dphi = TVector2::Phi_mpi_pi(-dphi);
	      if (reverse_dphi < -TMath::Pi()/2) reverse_dphi += 2*TMath::Pi();
	      h_deta_dphi->Fill(-deta, reverse_dphi, event_weight);
	      h_deta_dphi_cent[icent]->Fill(-deta, reverse_dphi, event_weight);
	      
	    }
	  }

	  //if(net_charge!=2 || net_baryon!=2)
	  //	  cout <<"net-charge = " << net_charge << ", net-baryon = " << net_baryon << endl;
  } // end event looping


  bool yieldClosure = true, probabilityClosure = true, multiplicityClosure = true;
  for (int b = 0; b <= nBins+1; ++b) {
    const double denominator = h_Egamma_generator_all.GetBinContent(b);
    double sumG = 0, sumU = 0, sumP = 0;
    for (int c = 0; c < 3; ++c) {
      double numerator = generatorClass[c]->GetBinContent(b);
      sumG += numerator;
      sumU += upcClass[c]->GetBinContent(b);
      double value = 0, variance = 0;
      if (denominator != 0) {
        value = numerator/denominator;
        if (unitGeneratorWeights) variance = value*(1-value)/denominator;
        else {
          // Subset covariance: numerator is part of denominator.
          double s2c = pow(generatorClass[c]->GetBinError(b), 2);
          double s2all = pow(h_Egamma_generator_all.GetBinError(b), 2);
          variance = ((1-2*value)*s2c + value*value*s2all)/(denominator*denominator);
        }
      } // Empty (or exactly cancelled) denominator: probability and error = 0.
      probability[c]->SetBinContent(b, value);
      probability[c]->SetBinError(b, sqrt(TMath::Max(0.0, variance)));
      sumP += value;
    }
    const double inclusiveUPC = h_Egamma_upc_all.GetBinContent(b);
    if (fabs(sumG-denominator) > 1e-9*(1+fabs(denominator)) ||
        fabs(sumU-inclusiveUPC) > 1e-9*(1+fabs(inclusiveUPC))) {
      yieldClosure = false; cerr << "Yield closure failed in energy bin " << b << endl;
    }
    if (denominator != 0 && fabs(sumP-1) > 1e-9) probabilityClosure = false;
  }
  // Include energy under/overflow: original logarithmic edges start at 10 MeV.
  TH1D *projectionGenerator = h2_Egamma_Ntof_generator.ProjectionY(
      "noonProjectionCheckGenerator", 0, nBins+1, "e");
  TH1D *projectionUpc = h2_Egamma_Ntof_upc.ProjectionY(
      "noonProjectionCheckUpc", 0, nBins+1, "e");
  for (int b = 0; b <= h_Ntof_generator.GetNbinsX()+1; ++b) {
    const double generatorValue = projectionGenerator->GetBinContent(b);
    const double upcValue = projectionUpc->GetBinContent(b);
    if (fabs(generatorValue-h_Ntof_generator.GetBinContent(b)) >
            1e-9*(1+fabs(generatorValue)) ||
        fabs(upcValue-h_Ntof_upc.GetBinContent(b)) >
            1e-9*(1+fabs(upcValue)) ||
        fabs(upcValue-h_npart_tof_matched->GetBinContent(b)) >
            1e-9*(1+fabs(upcValue))) multiplicityClosure = false;
  }
  delete projectionGenerator;
  delete projectionUpc;
  cout << "Total input events: " << nevents
       << "\nEvents with unresolved photons: " << hPhotonUnresolved.GetEntries()
       << "\nEvents with unresolved incoming protons: " << hProtonUnresolved.GetEntries()
       << "\nAccepted events: " << noonEvents.GetEntries() << endl;
  for (int c = 0; c < 3; ++c)
    cout << "Class " << c << ": events=" << classCount[c] << ", totalWeight sum=" << classWeight[c] << endl;
  cout << "Invalid/unset generator weights defaulted to 1: " << invalidGeneratorWeights << endl;
  if (flux.weights)
    cout << "nOOn flux fraction below PYTHIA coverage threshold " << coverageThresholdGeV
         << " GeV = " << fluxFractionBelowCoverage << " (over supplied flux-density support)" << endl;
  cout << "Class energy yield closure (generator and UPC, including flow bins): " << (yieldClosure ? "PASS" : "FAIL")
       << "\nProbability sum closure (nonzero denominator bins): " << (probabilityClosure ? "PASS" : "FAIL")
       << "\nGenerator/UPC TOF projection closure: "
       << (multiplicityClosure ? "PASS" : "FAIL") << endl;
  fout->cd();

  noonEvents.Write();
  h2_Egamma_Ntof_generator.Write(); h2_Egamma_Ntof_upc.Write();
  h2_Egamma_Nraw_generator.Write(); h2_Egamma_Nraw_upc.Write();
  h2_Egamma_Naccept_generator.Write(); h2_Egamma_Naccept_upc.Write();
  h2_Egamma_Ntpc_generator.Write(); h2_Egamma_Ntpc_upc.Write();
  h_Egamma_generator_all.Write(); h_Egamma_upc_all.Write();
  for (int c = 0; c < 3; ++c) {
    generatorClass[c]->Write(); upcClass[c]->Write(); probability[c]->Write();
  }
  for (int stage = 0; stage < 4; ++stage) {
    h_Egamma_stage_generator[stage]->Write();
    h_Egamma_stage_upc[stage]->Write();
  }
  p_meanNtof_vs_Egamma.Write();
  h_Nraw_generator.Write(); h_Nraw_upc.Write();
  h_Naccept_generator.Write(); h_Naccept_upc.Write();
  h_Ntpc_generator.Write(); h_Ntpc_upc.Write();
  h_Ntof_generator.Write(); h_Ntof_upc.Write();
  h_WgammaP_all.Write(); h_WgammaP_class0.Write(); h_WgammaP_class1.Write();
  hProtonUnresolved.Write();
  TNamed("noonPhotonSelection", "0 < q.E() < 5 GeV; egammaMeV=1000*egammaGeV; 0 <= -q.M2() < 0.01 GeV^2; W_gamma_p=(q+pIncoming).M() [GeV]").Write();
  TNamed("noonBranchUnits", "eventIndex: zero-based input chain entry; egammaGeV,wGammaPGeV: GeV; egammaMeV: MeV; q2GeV2: GeV^2; remaining branches: dimensionless counts, labels, or weights").Write();
  TNamed("noonTrackSelection", "status=1 charged PDG particle; |eta|<1.0; 0.2<pT<2.0 GeV/c; TPC-reconstructed tracks are analysis/pair tracks; TOF defines event activity only").Write();
  TNamed("noonDetectorResponse", Form("Sequential independent Bernoulli response after kinematic acceptance; tpcEfficiency=%.17g; tpcRandomSeed=%u; tofEfficiency=%.17g conditional on TPC; tofRandomSeed=%u; TRandom3", tpcEfficiency, tpcRandomSeed, tofEfficiency, tofRandomSeed)).Write();
  TNamed("noonMultiplicityClasses", "class0: Ntof<=5; class1: 5<Ntof<=15; class2: Ntof>15").Write();
  TNamed("noonPhotonDirection", "+1: q.Pz()>0; -1: q.Pz()<0; 0: q.Pz()==0; generator laboratory coordinate system").Write();
  TNamed("noonWeightConvention", "generatorWeight: finite scalar input weight, else 1; totalWeight=generatorWeight*noonFluxWeight; generator-labelled QA uses generatorWeight; UPC-labelled and legacy physics outputs use totalWeight exactly once").Write();
  TNamed("noonZeroGeneratorWeightPolicy", zeroGeneratorWeightIsUnset ? "zero scalar STAR generator weight treated as unset: use 1 (unweighted sample)" : "zero scalar generator weights retained").Write();
  TNamed("noonProbabilityConvention", "Zero denominator: probability=error=0 (undefined); unit generator weights: binomial errors; otherwise subset covariance with sum of squared weights; flow bins included").Write();
  TNamed("noonAfterburnerInterface", "For each noonEvents entry with multClass 0 or 1, call NeutronGenerator::GenerateEvent(egammaGeV) separately by class; no forward neutrons generated by this analysis; preserve class2 for closure").Write();
  if (useTextWeights) {
    TNamed("noonFluxWeightApplied", Form("yes; exact-bin text file=%s; x in MeV; zero table entries give zero UPC weight; no interpolation", textWeightPath.Data())).Write();
    hNoonTextWeight.Write("noonFluxWeightInput");
  } else if (flux.weights) {
    TNamed("noonFluxWeightApplied", Form("yes; file=%s; weight=%s; physical flux density=%s; x in MeV, density in 1/MeV; coverage threshold=%.17g GeV; integrals over supplied support", argv[4], argv[5], argv[6], coverageThresholdGeV)).Write();
    flux.weights->Write("noonFluxWeightInput");
    flux.density->Write("noonFluxDensityInput");
  } else {
    TNamed("noonFluxWeightApplied", "no; noonFluxWeight=1; unweighted by UPC flux").Write();
  }
  hAcceptedEventWeight.Write();
  hAcceptedEvents.Write();
  hPhotonUnresolved.Write();
  for (int species = 0; species < 7; ++species) {
    hEtaSpecies[species]->Write();
    if (species < 6) hShiftedY[species]->Write();
  }
  gWeights->Write();
  fullq2hist.Write();
  q2hist.Write();
  egamma.Write();
  ypionplus.Write();
  ykaonplus.Write();
  yprotonplus.Write();
  ypionmins.Write();
  ykaonmins.Write();
  yprotonmins.Write();
  yneutron.Write();

  ptpionplus.Write();
  ptkaonplus.Write();
  ptprotonplus.Write();

  ptpionmins.Write();
  ptkaonmins.Write();
  ptprotonmins.Write();

  ptneutron.Write();


  h_dNdEta->Write();
  h_dNdEta2->Write();
  h_pt_all->Write();

  h_npart_raw->Write();
  h_npart_tof_matched->Write();

  h_deta_dphi->Write();
  for (int i = 0; i < nCentBins; ++i) {
      h_dNdEta_cent[i]->Write();
      h_dNdEta2_cent[i]->Write();
      h_deta_dphi_cent[i]->Write();
  }

  TDirectory* subD = fout->mkdir(simType.Data());
  subD->cd();
  hEventCounter->Write();
  for(int i=0; i<10; i++)
  {
	  hPartEta[i]->Write();
  }


  // ROOT owns directory-attached objects, including stack histograms. Detach
  // before closing, then validate a fresh read of the actual file on disk.
  vector<TObject*> outputObjects;
  TIter outputIterator(fout->GetList());
  TObject *outputObject;
  while ((outputObject = outputIterator())) outputObjects.push_back(outputObject);
  for (size_t j = 0; j < outputObjects.size(); ++j) {
    TH1 *h = dynamic_cast<TH1*>(outputObjects[j]);
    TTree *tree = dynamic_cast<TTree*>(outputObjects[j]);
    if (h) h->SetDirectory(0);
    if (tree) tree->SetDirectory(0);
  }
  fout->Close();
  TFile validationFile(outFileName.Data(), "READ");
  gROOT->cd();
  const bool finiteClosure = !validationFile.IsZombie() && finiteOutput(&validationFile);
  validationFile.Close();
  cout << "Finite output branches and histogram contents/errors: " << (finiteClosure ? "PASS" : "FAIL") << endl;

  //delete gWeights;

  cout<<"end of program"<<endl;
  return (yieldClosure && probabilityClosure && multiplicityClosure && finiteClosure) ? 0 : 2;
  exit(0);
}

//============================================
bool isPrimaryParticle(const int id)
{
	return (fabs(id)==211 || fabs(id)==321 || fabs(id)==2212 || fabs(id)==11 || fabs(id)==13);
}

bool isChargedParticle(const int id)
{
	return (abs(id)==11 || abs(id)==13 || abs(id)==211 ||
	        abs(id)==321 || abs(id)==2212);
}
