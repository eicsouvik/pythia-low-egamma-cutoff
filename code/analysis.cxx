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
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
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

// define histogramshPartEta
TH1F *hEventCounter;
TH1F *hPartEta[10];

int main(int argc, char **argv)
{
  if(argc < 4 ) {
    cerr << "Usage: " << argv[0] << " <simulation type> <file list> <output file>" << endl;
    return 1;
  }
  TString simType = argv[1];
  TString fileList = argv[2];
  TString outFileName = argv[3];
  
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

  TH1::SetDefaultSumw2();
  hEventCounter = new TH1F("hEventCounter", "Event statistics;sub-process", 110, 0, 110);

  const int nBins = 300;
  double logBins[nBins + 1];
  const double log_kMin = log10(10.0), log_kMax = log10(5e6);

  for (int i = 0; i <= nBins; ++i) {
         logBins[i]  = pow(10, log_kMin + (log_kMax - log_kMin) * i / nBins);
  }

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
  TH1D *h_dNdEta2 = new TH1D("h_dNdEta2", "Selected charged tracks, |#eta|<1, 0.2<p_{T}<2.0;#eta;Weighted yield", 201, -5.025, 5.025);
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

  TH1D *h_pt_all = new TH1D("h_pt_all", "p_{T} for all final state particles;p_{T} [GeV/c];Counts", 200, 0., 5.);

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
        Form("Selected charged tracks, |#eta|<1, 0.2<p_{T}<2.0 %s;#eta;Weighted yield", centBinTitles[i]),
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

  TGraph *gWeights = new TGraph();
  gWeights->SetName("gWeights");
  gWeights->SetTitle("Event weights vs E_{#gamma};E_{#gamma} [MeV];Weight");
  std::ifstream weightFile("reweighting_factors_ep2x100_test2.txt");
  std::string line;
  double energy_mev, weight;
  int point_index;
  int n_weights = 0;
  if (weightFile.is_open()) {
    std::getline(weightFile, line);
    while (std::getline(weightFile, line)) {
        std::stringstream ss(line);
        if (ss >> point_index >> energy_mev >> weight) {
            gWeights->SetPoint(n_weights, energy_mev, weight);
            ++n_weights;
        }
    }
    weightFile.close();
    std::cout << "Read " << n_weights << " points from reweighting_factors_ep2x100_test2.txt" << std::endl;
  } else {
      std::cerr << "Warning: Could not open reweighting_factors_ep2x100_test2.txt; weights will be 1." << std::endl;
      gWeights->SetPoint(0, 0.0, 1.0);
      gWeights->SetPoint(1, 5000.0, 1.0);
  }


  if (n_weights == 0) {
    std::cerr << "Warning: No valid weights found; weights will be 1." << std::endl;
    gWeights->SetPoint(0, 0.0, 1.0);
    gWeights->SetPoint(1, 5000.0, 1.0);
  }

  //------------------------------------------------
  // Loop the chain
  //------------------------------------------------
  int nevents = (int)chain->GetEntries();
  printf("\n----------------------------------\n");
  printf("Loop the chain: %d events\n",nevents);
  printf("----------------------------------\n\n");

  genevents *t = new genevents(chain);

  for(int i=0;i<nevents;i++)
  {
	  //if(i==2000) break;
	  if(i%1000==0) cout << "begin " << i << "th entry...." << endl;
	  t->GetEntry(i);

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


          TLorentzVector q;
          if (!exchangedPhoton(t, q)) {
            hPhotonUnresolved.Fill(0.5);
            continue;
          }
          const double E_gamma = q.E(); // beam/laboratory frame
          const double Q2_gamma = -q.M2(); // Q2fac is NOT photon virtuality
          const double event_weight = 1.0;//gWeights->Eval(E_gamma * 1e3); // Convert GeV to MeV for the loaded factors.
          fullq2hist.Fill(Q2_gamma, event_weight);
          
	  if (!(E_gamma > 0.0 && E_gamma < 5.0)) continue;
          if (!(Q2_gamma >= 0.0 && Q2_gamma < 0.01)) continue;
          
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

	  // const double tof_efficiency = 0.9;
	  const double track_eta_max = 1.0;
	  const double track_pt_min = 0.2;
	  const double track_pt_max = 2.0;
	  std::vector<TLorentzVector> selected_tracks;

	  // Apply acceptance to each eligible track; TOF efficiency rejection is disabled.
	  // A sampled count must not be used as a PYTHIA event-record boundary.
	  for (int j = 0; j < npart; ++j) {
	    if (t->mParticles_mStatus[j] != 1) continue;
	    if (!isChargedParticle(t->mParticles_mId[j])) continue;

	    TLorentzVector track;
	    track.SetPxPyPzE(t->mParticles_mPx[j], t->mParticles_mPy[j],
	                    t->mParticles_mPz[j], t->mParticles_mEnergy[j]);

	    if (fabs(track.Eta()) >= track_eta_max) continue;
	    if (track.Pt() <= track_pt_min || track.Pt() >= track_pt_max) continue;
	    // if (gRandom->Rndm() >= tof_efficiency) continue;
	    selected_tracks.push_back(track);
	  }

	  int npart_tof_matched = selected_tracks.size();
	  h_npart_tof_matched->Fill(npart_tof_matched, event_weight);

	  int icent = -1;
          if (npart_tof_matched <= 5) {
            icent = 0;
	  } else if (npart_tof_matched > 5 && npart_tof_matched <= 15) {
            icent = 1;
	  } else {
	    icent = 2;
	  }


          for (const auto& track : selected_tracks) {
            h_dNdEta2->Fill(track.Eta(), event_weight);
            h_dNdEta2_cent[icent]->Fill(track.Eta(), event_weight);
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
		  h_pt_all->Fill(pt, event_weight);


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
	  for (size_t j = 0; j < selected_tracks.size(); ++j) {
	    for (size_t k = j + 1; k < selected_tracks.size(); ++k) {
	      double deta = selected_tracks[k].Eta() - selected_tracks[j].Eta();
	      double dphi = TVector2::Phi_mpi_pi(
	          selected_tracks[k].Phi() - selected_tracks[j].Phi());

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


  TFile *fout = new TFile(outFileName.Data(),"RECREATE");
  fout->cd();

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


  fout->Close();
  //delete gWeights;

  cout<<"end of program"<<endl;
  return(0);
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
