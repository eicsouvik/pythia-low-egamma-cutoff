#include "NeutronGenerator.h"
#include "TRandom.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>

// Run from noon_original. maxEntries=-1 processes the complete tree.
// Energy cuts are optional, in GeV; the activity classes do NOT define disjoint k windows.
void run_pythia_noon(
 const char* input="../PYTHIA_root_files/egamma_less5gev_ep10x100_test17-18_updated_analysis_macro2.root",
 Long64_t maxEntries=-1,
 const char* outDir="pythia_noon_output",
 UInt_t seed=12345, double kMin=0., double kMax=5.,
 Long64_t firstEntry=0, Long64_t checkpointEvery=100000)
{
 try {
  if (!(kMax>kMin) || maxEntries==0 || firstEntry<0 || checkpointEvery<1) throw std::runtime_error("Invalid energy range, entry limit, offset or checkpoint interval");
  std::unique_ptr<TFile> in(TFile::Open(input));
  if(!in || in->IsZombie()) throw std::runtime_error("Cannot open input");
  auto tree=dynamic_cast<TTree*>(in->Get("noonEvents"));
  if(!tree) throw std::runtime_error("Missing noonEvents");
  if(firstEntry>=tree->GetEntries()) throw std::runtime_error("firstEntry is beyond the input tree");
  const char* required[]={"eventIndex","egammaGeV","nTofMatched","multClass","generatorWeight","noonFluxWeight","totalWeight","photonDirection"};
  for(auto b:required) if(!tree->GetBranch(b)) throw std::runtime_error(TString::Format("Missing branch %s",b).Data());
  if(gSystem->AccessPathName("Data/ENDF_Au197.txt")) throw std::runtime_error("Run from noon_original (Data/ENDF_Au197.txt missing)");
  gSystem->mkdir(outDir,kTRUE);
  TString output=TString::Format("%s/pythia_noon.root",outDir);
  if(!gSystem->AccessPathName(output)) throw std::runtime_error("Output already exists: choose a fresh outDir");
  TH1::AddDirectory(kFALSE);
  gRandom->SetSeed(seed);
  // No FinishProduction: it writes a fixed Output.root and unweighted internal QA.
  // Our weighted output is written explicitly below. Avoid legacy destructor ownership issues.
  auto gen=new NeutronGenerator;
  gen->SetBeamParameters(NeutronGenerator::kAu197,107.3);
  gen->SetDataPath("./Data/");
  gen->SetRunMode(NeutronGenerator::kInterface);
  gen->Initialize(); gen->Setup();

  TFile out(output,"RECREATE");
  if(out.IsZombie()) throw std::runtime_error("Cannot create output file");
  TH1D *energyIn[2],*energyOut[2],*neutrons[2],*plus[2],*minus[2],*ntof[2],*tags[2];
  TH2D *sides[2],*energyN[2];
  for(int c=0;c<2;c++) {
   energyIn[c]=new TH1D(Form("hEgammaInput_c%d",c),";E_{#gamma} [GeV];Weighted events",100,0,5);
   energyOut[c]=new TH1D(Form("hEgammaAfterburner_c%d",c),";E_{#gamma} [GeV];Weighted events",100,0,5);
   neutrons[c]=new TH1D(Form("hNeutronTotal_c%d",c),";Generated neutron count (both beams);Weighted events",101,-.5,100.5);
   plus[c]=new TH1D(Form("hNeutronPlus_c%d",c),";Generated neutrons, p_{z}>0;Weighted events",51,-.5,50.5);
   minus[c]=new TH1D(Form("hNeutronMinus_c%d",c),";Generated neutrons, p_{z}<0;Weighted events",51,-.5,50.5);
   ntof[c]=new TH1D(Form("hNtof_c%d",c),";N_{trk}^{TOF};Weighted events",17,-.5,16.5);
   tags[c]=new TH1D(Form("hNeutronTags_c%d",c),";Neutron tag;Weighted events",3,.5,3.5);
   tags[c]->GetXaxis()->SetBinLabel(1,"0n0n");tags[c]->GetXaxis()->SetBinLabel(2,"0nXn + Xn0n");tags[c]->GetXaxis()->SetBinLabel(3,"XnXn");
   sides[c]=new TH2D(Form("hNeutronSides_c%d",c),";N_{n}(p_{z}>0);N_{n}(p_{z}<0)",51,-.5,50.5,51,-.5,50.5);
   energyN[c]=new TH2D(Form("hEgammaNeutrons_c%d",c),";E_{#gamma} [GeV];Generated neutron count",100,0,5,101,-.5,100.5);
   for(auto h:{energyIn[c],energyOut[c],neutrons[c],plus[c],minus[c],ntof[c],tags[c]}) h->Sumw2();
   sides[c]->Sumw2();energyN[c]->Sumw2();
   for(auto h:{energyIn[c],energyOut[c],neutrons[c],plus[c],minus[c],ntof[c],tags[c]})h->SetStats(kFALSE);
  }
  TTree events("neutronEvents","PYTHIA event linked to nOOn electromagnetic breakup");
  events.SetDirectory(&out);
  events.SetAutoFlush(10000); // Fixed basket flush interval; independent of compression.
  events.SetAutoSave(0); // Explicit coordinated checkpoints below.
  Long64_t eventIndex; int cls,nTof,nPlus,nMinus,direction; double k,w,wgen,wflux;
  events.Branch("eventIndex",&eventIndex);events.Branch("multClass",&cls);
  events.Branch("egammaGeV",&k);events.Branch("nTofMatched",&nTof);
  events.Branch("photonDirection",&direction);events.Branch("nNeutronsPlus",&nPlus);events.Branch("nNeutronsMinus",&nMinus);
  events.Branch("generatorWeight",&wgen);events.Branch("noonFluxWeight",&wflux);events.Branch("totalWeight",&w);
  TTreeReader reader(tree);
  reader.SetEntriesRange(firstEntry,tree->GetEntries());
  TTreeReaderValue<Long64_t> idx(reader,"eventIndex");
  TTreeReaderValue<Double_t> ek(reader,"egammaGeV"),wg(reader,"generatorWeight"),wf(reader,"noonFluxWeight"),wt(reader,"totalWeight");
  TTreeReaderValue<Int_t> mc(reader,"multClass"),nm(reader,"nTofMatched"),pd(reader,"photonDirection");
  Long64_t scanned=0,counts[2]={0,0},rejected=0; double sums[2]={0,0};
  auto checkpoint=[&]() {
   out.cd();
   if(events.AutoSave("SaveSelf")<0) throw std::runtime_error("Tree checkpoint failed");
   for(int c=0;c<2;c++) {
    for(auto h:{energyIn[c],energyOut[c],neutrons[c],plus[c],minus[c],ntof[c],tags[c]})h->Write("",TObject::kOverwrite);
    sides[c]->Write("",TObject::kOverwrite);energyN[c]->Write("",TObject::kOverwrite);
    TParameter<Long64_t>(Form("generatedEvents_c%d",c),counts[c]).Write("",TObject::kOverwrite);
    TParameter<double>(Form("sumWeights_c%d",c),sums[c]).Write("",TObject::kOverwrite);
   }
   TParameter<Long64_t>("firstEntry",firstEntry).Write("",TObject::kOverwrite);
   TParameter<Long64_t>("nextInputEntry",firstEntry+scanned).Write("",TObject::kOverwrite);
   TParameter<Long64_t>("scannedEntries",scanned).Write("",TObject::kOverwrite);
   TParameter<Long64_t>("invalidEntries",rejected).Write("",TObject::kOverwrite);
   TParameter<Long64_t>("checkpointEvery",checkpointEvery).Write("",TObject::kOverwrite);
   TParameter<Long64_t>("randomSeed",seed).Write("",TObject::kOverwrite);
   TParameter<double>("kMinGeV",kMin).Write("",TObject::kOverwrite);
   TParameter<double>("kMaxGeV",kMax).Write("",TObject::kOverwrite);
   TNamed("inputFile",input).Write("",TObject::kOverwrite);
   TNamed("classDefinition","class0: Ntof<=5; class1: 5<Ntof<=15; overlapping conditional photon spectra").Write("",TObject::kOverwrite);
   for(auto name:{"noonFluxWeightApplied","noonPhotonSelection","noonTrackSelection","noonMultiplicityClasses","noonTofResponse"}) if(auto o=in->Get(name))o->Write(name,TObject::kOverwrite);
   TParameter<int>("runComplete",0).Write("",TObject::kOverwrite);
   out.SaveSelf(kTRUE);out.Flush();
   if(out.TestBit(TFile::kWriteError)) throw std::runtime_error("Output checkpoint write error");
   std::cout<<"Checkpoint saved: "<<scanned<<" scanned; next input entry "<<firstEntry+scanned<<std::endl;
  };
  checkpoint();
  while((maxEntries<0 || scanned<maxEntries) && reader.Next()) {
   if(scanned>0 && scanned%checkpointEvery==0) checkpoint();
   ++scanned; cls=*mc; if(cls!=0 && cls!=1) continue;
   k=*ek;w=*wt;wgen=*wg;wflux=*wf;nTof=*nm;eventIndex=*idx;direction=*pd;
   if((cls==0 && !(nTof>=0 && nTof<=5)) || (cls==1 && !(nTof>5 && nTof<=15))) throw std::runtime_error("Multiplicity class inconsistent with nTofMatched");
   if(!std::isfinite(k)||!std::isfinite(w)||!std::isfinite(wgen)||!std::isfinite(wflux)) {++rejected;continue;}
   if(k<=kMin || k>=kMax) continue;
   if(std::abs(w-wgen*wflux)>1.e-9*(1.+std::abs(w))) throw std::runtime_error("Inconsistent weight product");
   energyIn[cls]->Fill(k,w);ntof[cls]->Fill(nTof,w);
   // Original CreateNeutrons repeatedly calls ProjectionY("hENDF_1D").
   // Enable directory registration during generation so ROOT reuses this
   // projection instead of leaking one histogram per beam side per event.
   out.cd();TH1::AddDirectory(kTRUE);
   gen->GenerateEvent(k);
   TH1::AddDirectory(kFALSE);
   nPlus=0;nMinus=0;
   auto particles=gen->ImportParticles();
   for(int j=0;j<particles->GetEntriesFast();j++) {
    auto p=static_cast<TParticle*>(particles->At(j));
    if(p->GetPdgCode()!=2112) continue;
    if(p->Pz()>=0) ++nPlus;else ++nMinus;
   }
   energyOut[cls]->Fill(k,w);neutrons[cls]->Fill(nPlus+nMinus,w);
   plus[cls]->Fill(nPlus,w);minus[cls]->Fill(nMinus,w);
   sides[cls]->Fill(nPlus,nMinus,w);energyN[cls]->Fill(k,nPlus+nMinus,w);
   tags[cls]->Fill((nPlus==0 && nMinus==0)?1:((nPlus==0 || nMinus==0)?2:3),w);
   if(events.Fill()<0) throw std::runtime_error("Failed to write neutron event");
   ++counts[cls];sums[cls]+=w;
   gen->FinishEvent();
   if(scanned%10000==0) std::cout<<"Scanned "<<scanned<<" entries; generated "<<counts[0]<<", "<<counts[1]<<std::endl;
  }
  if(reader.GetEntryStatus()!=TTreeReader::kEntryBeyondEnd && firstEntry+scanned<tree->GetEntries() && (maxEntries<0 || scanned<maxEntries)) throw std::runtime_error("Tree reader stopped unexpectedly");
  checkpoint();
  out.cd();events.Write("",TObject::kOverwrite);
  for(int c=0;c<2;c++) {
   for(auto h:{energyIn[c],energyOut[c],neutrons[c],plus[c],minus[c],ntof[c],tags[c]})h->Write();
   sides[c]->Write();energyN[c]->Write();
   TParameter<Long64_t>(Form("generatedEvents_c%d",c),counts[c]).Write();
   TParameter<double>(Form("sumWeights_c%d",c),sums[c]).Write();
  }
  TParameter<Long64_t>("scannedEntries",scanned).Write();TParameter<Long64_t>("invalidEntries",rejected).Write();
  TParameter<Long64_t>("randomSeed",seed).Write();
  TParameter<double>("kMinGeV",kMin).Write();TParameter<double>("kMaxGeV",kMax).Write();
  TNamed("physicsScope","Au197 gamma=107.3; original nOOn channel conventions retained, including grouped channels; generated counts are model output, exclusive physical neutron bins require validation; electromagnetic breakup only").Write();
  TNamed("energyConvention","egammaGeV passed unchanged; input metadata describes generator laboratory energy. Matching its frame to the AuAu collider frame is an analysis assumption.").Write();
  TNamed("weightConvention","Input totalWeight used once; no extra UPC flux applied; no detector smearing; beam sides use pz sign, not STAR East/West").Write();
  for(auto name:{"noonFluxWeightApplied","noonPhotonSelection","noonTrackSelection","noonMultiplicityClasses","noonTofResponse"}) if(auto o=in->Get(name)) o->Write(name);

  gStyle->SetOptStat(0);
  auto annotation=[](int c) {TLatex t;t.SetNDC();t.SetTextSize(.035);t.DrawLatex(.15,.87,c==0?"N_{trk}^{TOF} #leq 5":"5 < N_{trk}^{TOF} #leq 15");t.SetTextSize(.028);t.DrawLatex(.15,.81,"PYTHIA #gamma+p + nOOn Au electromagnetic breakup");};
  auto save=[&](TCanvas& canvas,const char* name){out.cd();canvas.Write();canvas.SaveAs(Form("%s/%s.pdf",outDir,name));};
  TCanvas energies("energy_closure","Photon energy input and afterburner",1200,550);energies.Divide(2,1);
  for(int c=0;c<2;c++){energies.cd(c+1);energyIn[c]->SetLineColor(kBlue+1);energyIn[c]->Draw("HIST");energyOut[c]->SetMarkerStyle(24);energyOut[c]->SetMarkerColor(kRed+1);energyOut[c]->Draw("E SAME");TLegend l(.49,.66,.88,.77);l.SetBorderSize(0);l.AddEntry(energyIn[c],"PYTHIA input","l");l.AddEntry(energyOut[c],"Passed to nOOn","p");l.DrawClone();annotation(c);}save(energies,"qa_energy_closure");
  TCanvas result("neutron_distributions","Neutrons by activity",1200,550);result.Divide(2,1);
  for(int c=0;c<2;c++){result.cd(c+1);neutrons[c]->SetLineColor(c==0?kBlue+1:kRed+1);neutrons[c]->Draw("E");annotation(c);}save(result,"neutron_distributions");
  for(int c=0;c<2;c++){
   TCanvas single(Form("neutron_class%d",c),"Neutron distribution",800,600);
   neutrons[c]->Draw("E");annotation(c);save(single,Form("neutron_distribution_class%d",c));
  }
  TCanvas shapes("class_shape_comparison","Activity-conditioned photon and neutron spectra",1200,550);shapes.Divide(2,1);
  std::unique_ptr<TH1D> shape[4];
  for(int pad=0;pad<2;pad++){
   shapes.cd(pad+1);
   for(int c=0;c<2;c++){
    auto source=pad==0?energyIn[c]:neutrons[c];
    shape[2*pad+c].reset(static_cast<TH1D*>(source->Clone(Form("shape_%d_%d",pad,c))));
    auto h=shape[2*pad+c].get();double sum=h->Integral(0,h->GetNbinsX()+1);
    if(sum>0)h->Scale(1./sum);
    h->GetYaxis()->SetTitle("Unit-normalized weighted events");h->SetLineColor(c==0?kBlue+1:kRed+1);
   }
   auto a=shape[2*pad].get();auto b=shape[2*pad+1].get();
   a->SetMaximum(1.3*std::max(a->GetMaximum(),b->GetMaximum()));a->Draw("HIST");b->Draw("HIST SAME");
   TLegend legend(.49,.72,.88,.87);legend.SetBorderSize(0);legend.AddEntry(a,"N_{trk}^{TOF} #leq 5","l");legend.AddEntry(b,"5 < N_{trk}^{TOF} #leq 15","l");legend.DrawClone();
  }
  save(shapes,"qa_class_shape_comparison");
  TCanvas correlations("energy_neutron_QA","Energy versus neutrons",1200,550);correlations.Divide(2,1);
  for(int c=0;c<2;c++){correlations.cd(c+1);energyN[c]->Draw("COLZ");annotation(c);}save(correlations,"qa_energy_neutrons");
  TCanvas beam("beam_neutron_QA","Beam neutron counts",1200,550);beam.Divide(2,1);
  for(int c=0;c<2;c++){beam.cd(c+1);sides[c]->Draw("COLZ");annotation(c);}save(beam,"qa_neutron_sides");
  TCanvas activity("activity_QA","PYTHIA activity",1200,550);activity.Divide(2,1);
  for(int c=0;c<2;c++){activity.cd(c+1);ntof[c]->Draw("E");annotation(c);}save(activity,"qa_pythia_activity");
  TCanvas fractions("tag_QA","Neutron tags",1200,550);fractions.Divide(2,1);
  for(int c=0;c<2;c++){fractions.cd(c+1);tags[c]->Draw("E");annotation(c);}save(fractions,"qa_neutron_tags");
  TParameter<int>("runComplete",1).Write("",TObject::kOverwrite);
  out.Close();
  std::cout<<"Saved "<<output<<"; class counts "<<counts[0]<<", "<<counts[1]<<std::endl;
 } catch(const std::exception& e) {std::cerr<<"run_pythia_noon failed: "<<e.what()<<std::endl;throw;}
}
