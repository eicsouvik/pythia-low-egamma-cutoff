#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TNamed.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace noon_plot {
const int colors[2]={kAzure+2,kOrange+7};
const char* classes[2]={"N_{trk}^{TOF} #leq 5","5 < N_{trk}^{TOF} #leq 15"};
void text(double x,double y,const char* s,double size=.037,int color=kBlack){TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(size);t.SetTextColor(color);t.DrawLatex(x,y,s);}
void pad(bool heat=false){gPad->SetLeftMargin(.14);gPad->SetRightMargin(heat?.18:.05);gPad->SetTopMargin(.24);gPad->SetBottomMargin(.15);gPad->SetTicks(1,1);}
void heading(const char* title,int cls=-1,bool neutron=false){
 text(.14,.95,title,.036);text(.14,.895,cls<0?"PYTHIA #gamma+p + nOOn Au electromagnetic breakup":classes[cls],.033);
 text(.14,.84,neutron?"Generated counts; original channel mapping retained":"Run input weights; no additional UPC-flux reweighting",.028,kGray+2);
}
void style(TH1* h,const char* x,const char* y,int color){h->SetTitle("");h->SetStats(0);h->SetLineWidth(2);h->SetLineColor(color);h->SetMarkerColor(color);h->SetMarkerStyle(20);h->SetMarkerSize(.55);h->GetXaxis()->SetTitle(x);TString label(y);if(label.Contains("neutron-count"))label="Probability / neutron bin";if(label.Contains("Sigma"))label="Normalized photon density [GeV^{-1}]";h->GetYaxis()->SetTitle(label);for(auto a:{h->GetXaxis(),h->GetYaxis()}){a->SetLabelSize(.038);a->SetTitleSize(.040);}h->GetXaxis()->SetTitleOffset(1.15);h->GetYaxis()->SetTitleOffset(1.5);}
void legend(TH1* a,const char* labelA,TH1* b,const char* labelB){TLegend l(.43,.59,.93,.74);l.SetTextFont(42);l.SetTextSize(.030);l.SetFillStyle(0);l.SetBorderSize(0);l.AddEntry(a,labelA,"lp");if(b)l.AddEntry(b,labelB,"lp");l.DrawClone();}
void range(TH1* a,TH1* b=nullptr){a->SetMinimum(0);double mx=a->GetMaximum();if(b)mx=std::max(mx,b->GetMaximum());a->SetMaximum(mx>0?mx*1.35:1.);}

TString dated(const char* pattern){
 TString result(pattern);if(!result.Contains("mm_dd_yy"))return result;
 TDatime now;result.ReplaceAll("mm_dd_yy",Form("%02d_%02d_%02d",now.GetMonth(),now.GetDay(),now.GetYear()%100));return result;
}

// The n00n input is stored as a density versus E_gamma in MeV, while the
// afterburner/PYTHIA histograms use GeV.  Integrate the n00n density over each
// PYTHIA bin and return an average density per GeV on the PYTHIA axis.
TH1D* mapFluxToPythiaAxis(const TH1* fluxMeV,const TH1* pythia,const char* name){
 const TAxis* axis=pythia->GetXaxis();TH1D* mapped=nullptr;
 if(axis->GetXbins()->GetSize()>0)mapped=new TH1D(name,"",axis->GetNbins(),axis->GetXbins()->GetArray());
 else mapped=new TH1D(name,"",axis->GetNbins(),axis->GetXmin(),axis->GetXmax());
 mapped->SetDirectory(nullptr);mapped->Sumw2();
 for(int b=1;b<=mapped->GetNbinsX();++b){
  const double lowGeV=mapped->GetXaxis()->GetBinLowEdge(b),highGeV=mapped->GetXaxis()->GetBinUpEdge(b);
  const double lowMeV=1000.*lowGeV,highMeV=1000.*highGeV;double integral=0,error2=0;
  for(int f=1;f<=fluxMeV->GetNbinsX();++f){
   const double flo=fluxMeV->GetXaxis()->GetBinLowEdge(f),fhi=fluxMeV->GetXaxis()->GetBinUpEdge(f);
   const double overlap=std::max(0.,std::min(highMeV,fhi)-std::max(lowMeV,flo));if(overlap<=0)continue;
   integral+=fluxMeV->GetBinContent(f)*overlap;const double term=fluxMeV->GetBinError(f)*overlap;error2+=term*term;
  }
  const double widthGeV=highGeV-lowGeV;mapped->SetBinContent(b,widthGeV>0?integral/widthGeV:0);mapped->SetBinError(b,widthGeV>0?std::sqrt(error2)/widthGeV:0);
 }
 return mapped;
}

void normalizeCountsDensity(TH1* h,int first,int last){
 for(int b=0;b<=h->GetNbinsX()+1;++b)if(b<first||b>last){h->SetBinContent(b,0);h->SetBinError(b,0);}
 const double sum=h->Integral(first,last);if(sum<=0)throw std::runtime_error("Cannot normalize empty PYTHIA spectrum");h->Scale(1./sum,"width");
}

void normalizeStoredDensity(TH1* h,int first,int last){
 for(int b=0;b<=h->GetNbinsX()+1;++b)if(b<first||b>last){h->SetBinContent(b,0);h->SetBinError(b,0);}
 const double sum=h->Integral(first,last,"width");if(sum<=0)throw std::runtime_error("Cannot normalize empty n00n spectrum");h->Scale(1./sum);
}
}

// Run: root -l -b -q 'plot_pythia_noon_qa.C+("pythia_noon.root","essential_QA_plots_mm_dd_yy","../n00n_flux_from_Ashik_da/noon_flux.root")'
// Uses saved histograms only; never reruns or changes neutron generation.
void plot_pythia_noon_qa(const char* input="pythia_noon.root",const char* output="essential_QA_plots_mm_dd_yy",const char* fluxInput="../n00n_flux_from_Ashik_da/noon_flux.root"){
 using namespace noon_plot;
 std::unique_ptr<TFile> file(TFile::Open(input));
 if(!file||file->IsZombie())throw std::runtime_error("Cannot open input ROOT file");
 std::unique_ptr<TFile> fluxFile(TFile::Open(fluxInput));
 if(!fluxFile||fluxFile->IsZombie())throw std::runtime_error("Cannot open n00n flux input ROOT file");
 TString outputPath=dated(output);TString rootPath=TString::Format("%s/qa_plots.root",outputPath.Data());
 if(!gSystem->AccessPathName(rootPath))throw std::runtime_error("Output exists; choose a fresh subdirectory");
 std::vector<std::unique_ptr<TH1>> owned;
 auto get=[&](const char* key)->TH1*{auto src=dynamic_cast<TH1*>(file->Get(key));if(!src)throw std::runtime_error(TString::Format("Missing histogram %s",key).Data());auto h=static_cast<TH1*>(src->Clone(TString::Format("plot_%s",key)));h->SetDirectory(nullptr);owned.emplace_back(h);return h;};
 auto copy=[&](TH1* src,const char* name)->TH1*{auto h=static_cast<TH1*>(src->Clone(name));h->SetDirectory(nullptr);owned.emplace_back(h);return h;};
 auto probability=[&](TH1* src,const char* name,bool density=false)->TH1*{auto h=copy(src,name);double sum=h->Integral(0,h->GetNbinsX()+1);for(int b=0;b<=h->GetNbinsX()+1;b++)if(h->GetBinContent(b)<0)throw std::runtime_error("Negative weights: probability/log plots need a different prescription");if(sum<=0)throw std::runtime_error("Cannot normalize empty histogram");h->Scale(1./sum);if(density)h->Scale(1.,"width");return h;};
 TH1 *energy[2],*passed[2],*total[2],*plus[2],*minus[2],*activity[2],*tags[2];TH2 *joint[2],*sides[2];
 for(int c=0;c<2;c++){energy[c]=get(Form("hEgammaInput_c%d",c));passed[c]=get(Form("hEgammaAfterburner_c%d",c));total[c]=get(Form("hNeutronTotal_c%d",c));plus[c]=get(Form("hNeutronPlus_c%d",c));minus[c]=get(Form("hNeutronMinus_c%d",c));activity[c]=get(Form("hNtof_c%d",c));tags[c]=get(Form("hNeutronTags_c%d",c));joint[c]=dynamic_cast<TH2*>(get(Form("hEgammaNeutrons_c%d",c)));sides[c]=dynamic_cast<TH2*>(get(Form("hNeutronSides_c%d",c)));if(!joint[c]||!sides[c])throw std::runtime_error("Invalid map types");}
 gSystem->mkdir(outputPath,kTRUE);TFile out(rootPath,"RECREATE");if(out.IsZombie())throw std::runtime_error("Cannot create output");
 gStyle->SetOptStat(0);gStyle->SetOptTitle(0);gStyle->SetPalette(kViridis);gStyle->SetNumberContours(60);
 int figures=0;
 auto save=[&](TCanvas& c,const char* name){out.cd();c.Write(name);c.SaveAs(Form("%s/%s.pdf",outputPath.Data(),name));c.SaveAs(Form("%s/%s.png",outputPath.Data(),name));++figures;};
 // Closure ratios compare the SAME events and weights: they are correlated.
 TCanvas closure("closure","Photon energy transfer QA",1400,1000);closure.Divide(2,2);
 for(int c=0;c<2;c++){
  closure.cd(c+1);pad();style(energy[c],"E_{#gamma} [GeV]","Weighted events / 0.05 GeV",colors[0]);style(passed[c],"E_{#gamma} [GeV]","Weighted events / 0.05 GeV",colors[1]);passed[c]->SetMarkerStyle(24);range(energy[c],passed[c]);energy[c]->Draw("HIST");passed[c]->Draw("E1 SAME");legend(energy[c],"PYTHIA tree input",passed[c],"Energy passed to nOOn");heading("Primary photon energy: interface closure",c);
  double maxDiff=0;for(int b=0;b<=energy[c]->GetNbinsX()+1;b++)maxDiff=std::max(maxDiff,std::abs(energy[c]->GetBinContent(b)-passed[c]->GetBinContent(b)));text(.18,.55,Form("Maximum bin difference = %.3g",maxDiff),.032);
  closure.cd(c+3);pad();auto ratio=copy(passed[c],Form("energyClosureRatio_c%d",c));ratio->Reset();
  for(int b=1;b<=ratio->GetNbinsX();b++)if(energy[c]->GetBinContent(b)!=0){ratio->SetBinContent(b,passed[c]->GetBinContent(b)/energy[c]->GetBinContent(b));ratio->SetBinError(b,0);}
  style(ratio,"E_{#gamma} [GeV]","nOOn input / PYTHIA input",colors[c]);ratio->SetMinimum(.98);ratio->SetMaximum(1.02);ratio->Draw("AXIS");TLine line(0,1,5,1);line.SetLineStyle(2);line.DrawClone();
  // Draw only occupied denominator bins; undefined ratios are omitted.
  for(int b=1;b<=ratio->GetNbinsX();b++)if(energy[c]->GetBinContent(b)!=0){TMarker marker(ratio->GetBinCenter(b),ratio->GetBinContent(b),20);marker.SetMarkerColor(colors[c]);marker.SetMarkerSize(.45);marker.DrawClone();}
  heading("Correlated closure ratio",c);text(.16,.68,"Zero-denominator bins omitted; no independent-sample errors",.028);
 }save(closure,"qa_photon_energy_closure");
 // Separate, consistently normalized class distributions.
 for(int c=0;c<2;c++){
  TCanvas panel(Form("class%d",c),"Activity-conditioned neutron output",1400,1000);panel.Divide(2,2);
  panel.cd(1);pad();auto e=probability(energy[c],Form("photonDensity_c%d",c),true);style(e,"E_{#gamma} [GeV]","(1 / #Sigmaw) d#Sigmaw / dE_{#gamma} [GeV^{-1}]",colors[c]);range(e);e->Draw("E1");heading("PYTHIA photon spectrum for this activity class",c);legend(e,Form("PYTHIA; #LT E_{#gamma} #GT = %.3f GeV",energy[c]->GetMean()),nullptr,"");
  panel.cd(2);pad();auto n=probability(total[c],Form("neutronProbability_c%d",c));style(n,"N_{n}^{+} + N_{n}^{-}","Weighted probability / neutron-count bin",colors[c]);range(n);n->Draw("E1");heading("nOOn generated total neutron distribution",c,true);legend(n,"nOOn (includes zero-neutron events)",nullptr,"");
  panel.cd(3);pad();auto zoom=copy(n,Form("neutronLowCount_c%d",c));zoom->GetXaxis()->SetRangeUser(-.5,15.5);range(zoom);zoom->Draw("E1");heading("Low-count detail: original normalization",c,true);text(.17,.67,"Zoom only; events outside this range are not renormalized",.028);
  panel.cd(4);pad();auto p=probability(plus[c],Form("beamPlusProbability_c%d",c));auto m=probability(minus[c],Form("beamMinusProbability_c%d",c));style(p,"Generated neutron count on one beam side","Weighted probability / neutron-count bin",colors[0]);style(m,"Generated neutron count on one beam side","Weighted probability / neutron-count bin",colors[1]);m->SetMarkerStyle(24);range(p,m);p->Draw("E1");m->Draw("E1 SAME");legend(p,"nOOn: p_{z} > 0",m,"nOOn: p_{z} < 0");heading("Beam-side QA (not STAR East/West labels)",c,true);
  save(panel,Form("neutrons_and_photon_spectrum_class%d",c));
  TCanvas single(Form("total_class%d",c),"Total neutron spectrum",950,750);pad();n->Draw("E1");heading("Total neutron distribution",c,true);legend(n,classes[c],nullptr,"");save(single,Form("neutron_distribution_class%d",c));
 }
 TCanvas comparison("comparison","Conditional class spectra",1400,600);comparison.Divide(2,1);
 for(int p=0;p<2;p++){comparison.cd(p+1);pad();TH1* h[2];for(int c=0;c<2;c++){h[c]=probability(p==0?energy[c]:total[c],Form("compare_%d_c%d",p,c),p==0);style(h[c],p==0?"E_{#gamma} [GeV]":"Generated total neutron count",p==0?"Normalized photon density [GeV^{-1}]":"Weighted probability / neutron-count bin",colors[c]);}range(h[0],h[1]);h[0]->Draw("HIST");h[1]->Draw("HIST SAME");legend(h[0],classes[0],h[1],classes[1]);heading(p==0?"Overlapping activity-conditioned photon spectra":"Activity-conditioned neutron shapes",-1,p==1);}save(comparison,"qa_class_shape_comparison");
 TCanvas activityCanvas("activity","PYTHIA activity and nOOn tags",1400,1000);activityCanvas.Divide(2,2);
 for(int c=0;c<2;c++){activityCanvas.cd(c+1);pad();auto a=probability(activity[c],Form("activityProbability_c%d",c));style(a,"N_{trk}^{TOF}","Probability within activity class",colors[c]);a->GetXaxis()->SetRangeUser(c==0?-.5:5.5,c==0?5.5:15.5);range(a);a->Draw("E1");heading("PYTHIA TOF activity selection",c);legend(a,Form("#LT N_{trk}^{TOF} #GT = %.3f",activity[c]->GetMean()),nullptr,"");activityCanvas.cd(c+3);pad();auto t=probability(tags[c],Form("tagProbability_c%d",c));style(t,"Generated neutron tag","Fraction of weighted events",colors[c]);t->SetMinimum(0);t->SetMaximum(1);t->Draw("E1");heading("nOOn neutron-tag fractions (provisional)",c,true);for(int b=1;b<=3;b++){TLatex label;label.SetTextAlign(21);label.SetTextSize(.035);label.DrawLatex(b,t->GetBinContent(b)+.035,Form("%.2f%%",100*t->GetBinContent(b)));}}
 save(activityCanvas,"qa_activity_and_neutron_tags");
 // Maps with equal scales, and a conditional alternative removing the energy prior.
 for(int mode=0;mode<3;mode++){
  TCanvas maps(Form("maps%d",mode),"Neutron correlation QA",1400,650);maps.Divide(2,1);TH2* h[2];double maximum=0;
  for(int c=0;c<2;c++){
   h[c]=static_cast<TH2*>(copy(mode==2?sides[c]:joint[c],Form("map_%d_c%d",mode,c)));
   if(mode==1){for(int x=1;x<=h[c]->GetNbinsX();x++){double sum=0;for(int y=0;y<=h[c]->GetNbinsY()+1;y++)sum+=h[c]->GetBinContent(x,y);for(int y=0;y<=h[c]->GetNbinsY()+1;y++){h[c]->SetBinContent(x,y,sum>0?h[c]->GetBinContent(x,y)/sum:0);h[c]->SetBinError(x,y,0);}}}
   else {double sum=h[c]->Integral(0,h[c]->GetNbinsX()+1,0,h[c]->GetNbinsY()+1);if(sum<=0)throw std::runtime_error("Empty map");h[c]->Scale(1./sum);}
   maximum=std::max(maximum,h[c]->GetMaximum());
  }
  for(int c=0;c<2;c++){maps.cd(c+1);pad(true);gPad->SetLogz();style(h[c],mode==2?"N_{n}(p_{z}>0)":"E_{#gamma} [GeV]",mode==2?"N_{n}(p_{z}<0)":"Generated total neutron count",colors[c]);h[c]->SetMinimum(1.e-6);h[c]->SetMaximum(std::max(maximum,1.e-5));h[c]->GetZaxis()->SetTitle(mode==1?"P(N_{n} | E_{#gamma} bin)":"Probability / 2D bin");h[c]->GetZaxis()->SetTitleSize(.033);h[c]->GetZaxis()->SetLabelSize(.029);h[c]->GetZaxis()->SetTitleOffset(1.4);h[c]->Draw("COLZ");heading(mode==0?"Joint photon energy and neutron count":mode==1?"Neutron probability at fixed photon-energy bin":"Two-beam neutron probability",c,true);text(.15,.79,"Shared log scale; blank cells may have zero entries",.026);}
  save(maps,mode==0?"qa_energy_neutrons_joint":mode==1?"qa_energy_neutrons_conditional":"qa_neutron_beam_correlation");
 }

 // Compact overview of the essential objects stored in pythia_noon.root.
 TCanvas essentials("essentials","Essential PYTHIA+n00n outputs",1500,1050);essentials.Divide(2,2);
 essentials.cd(1);pad();TH1* essentialEnergy[2];
 for(int c=0;c<2;c++){essentialEnergy[c]=probability(energy[c],Form("essentialEnergy_c%d",c),true);style(essentialEnergy[c],"E_{#gamma} [GeV]","Normalized density [GeV^{-1}]",colors[c]);}
 range(essentialEnergy[0],essentialEnergy[1]);essentialEnergy[0]->Draw("HIST");essentialEnergy[1]->Draw("HIST SAME");legend(essentialEnergy[0],classes[0],essentialEnergy[1],classes[1]);heading("Activity-conditioned PYTHIA photon spectra");
 essentials.cd(2);pad();TH1* essentialActivity[2];
 for(int c=0;c<2;c++){essentialActivity[c]=probability(activity[c],Form("essentialActivity_c%d",c));style(essentialActivity[c],"N_{trk}^{TOF}","Probability",colors[c]);}
 range(essentialActivity[0],essentialActivity[1]);essentialActivity[0]->Draw("E1");essentialActivity[1]->Draw("E1 SAME");legend(essentialActivity[0],classes[0],essentialActivity[1],classes[1]);heading("PYTHIA event-activity distributions");
 essentials.cd(3);pad();TH1* essentialNeutrons[2];
 for(int c=0;c<2;c++){essentialNeutrons[c]=probability(total[c],Form("essentialNeutrons_c%d",c));style(essentialNeutrons[c],"N_{n}^{+}+N_{n}^{-}","Probability / neutron bin",colors[c]);}
 range(essentialNeutrons[0],essentialNeutrons[1]);essentialNeutrons[0]->Draw("HIST");essentialNeutrons[1]->Draw("HIST SAME");legend(essentialNeutrons[0],classes[0],essentialNeutrons[1],classes[1]);heading("n00n total-neutron distributions",-1,true);
 essentials.cd(4);pad();TH1* essentialTags[2];
 for(int c=0;c<2;c++){essentialTags[c]=probability(tags[c],Form("essentialTags_c%d",c));style(essentialTags[c],"Generated neutron tag","Fraction of weighted events",colors[c]);}
 essentialTags[0]->SetMinimum(0);essentialTags[0]->SetMaximum(1);essentialTags[0]->Draw("E1");essentialTags[1]->Draw("E1 SAME");legend(essentialTags[0],classes[0],essentialTags[1],classes[1]);heading("n00n neutron-tag fractions",-1,true);
 save(essentials,"qa_essential_overview");

 // Compare the inclusive PYTHIA output spectrum with the independent n00n
 // photon-flux input.  PYTHIA is stored in GeV; n00n is stored in MeV.
 TH2* flux2D=dynamic_cast<TH2*>(fluxFile->Get("h2fluxNoon_wPnohad"));
 if(!flux2D)throw std::runtime_error("Missing h2fluxNoon_wPnohad in n00n flux input");
 TH1D* fluxProjection=flux2D->ProjectionX("n00nFluxProjectionMeV",1,flux2D->GetNbinsY(),"e");fluxProjection->SetDirectory(nullptr);owned.emplace_back(fluxProjection);
 TH1* pythiaInclusive=copy(energy[0],"pythiaInclusiveEgamma");pythiaInclusive->Add(energy[1]);
 TH1D* fluxOnPythia=mapFluxToPythiaAxis(fluxProjection,pythiaInclusive,"n00nFluxOnPythiaAxis");owned.emplace_back(fluxOnPythia);
 int firstCommon=-1,lastCommon=-1;
 for(int b=1;b<=pythiaInclusive->GetNbinsX();++b){if(pythiaInclusive->GetBinContent(b)>0&&fluxOnPythia->GetBinContent(b)>0){if(firstCommon<0)firstCommon=b;lastCommon=b;}}
 if(firstCommon<0||lastCommon<firstCommon)throw std::runtime_error("No common E_gamma support between n00n and PYTHIA");
 normalizeCountsDensity(pythiaInclusive,firstCommon,lastCommon);normalizeStoredDensity(fluxOnPythia,firstCommon,lastCommon);
 TH1* pythiaOverFlux=copy(pythiaInclusive,"pythiaOverN00nFlux");pythiaOverFlux->Reset("ICES");
 for(int b=firstCommon;b<=lastCommon;++b){const double denominator=fluxOnPythia->GetBinContent(b);if(denominator<=0)continue;pythiaOverFlux->SetBinContent(b,pythiaInclusive->GetBinContent(b)/denominator);pythiaOverFlux->SetBinError(b,0);}
 const double commonLow=pythiaInclusive->GetXaxis()->GetBinLowEdge(firstCommon),commonHigh=pythiaInclusive->GetXaxis()->GetBinUpEdge(lastCommon);
 TCanvas fluxComparison("fluxComparison","n00n flux versus PYTHIA",1400,650);fluxComparison.Divide(2,1);
 fluxComparison.cd(1);pad();gPad->SetLogy();style(pythiaInclusive,"E_{#gamma} [GeV]","Unit-normalized density [GeV^{-1}]",colors[0]);style(fluxOnPythia,"E_{#gamma} [GeV]","Unit-normalized density [GeV^{-1}]",colors[1]);fluxOnPythia->SetMarkerStyle(24);
 double positiveMinimum=1.e300,maximum=0;for(auto h:{pythiaInclusive,(TH1*)fluxOnPythia})for(int b=firstCommon;b<=lastCommon;++b){const double value=h->GetBinContent(b);if(value>0)positiveMinimum=std::min(positiveMinimum,value);maximum=std::max(maximum,value);}pythiaInclusive->GetXaxis()->SetRangeUser(commonLow,commonHigh);pythiaInclusive->SetMinimum(positiveMinimum<1.e299?positiveMinimum*.4:1.e-12);pythiaInclusive->SetMaximum(maximum>0?maximum*5:1);pythiaInclusive->Draw("HIST");fluxOnPythia->Draw("HIST SAME");legend(pythiaInclusive,"PYTHIA output (classes 0+1)",fluxOnPythia,"n00n photon-flux input");heading("Photon-energy shape comparison");text(.17,.78,Form("Common support: %.3f < E_{#gamma} < %.3f GeV",commonLow,commonHigh),.027);
 fluxComparison.cd(2);pad();gPad->SetLogy();style(pythiaOverFlux,"E_{#gamma} [GeV]","PYTHIA / n00n",kMagenta+2);pythiaOverFlux->GetXaxis()->SetRangeUser(commonLow,commonHigh);double ratioMin=1.e300,ratioMax=0;for(int b=firstCommon;b<=lastCommon;++b){const double value=pythiaOverFlux->GetBinContent(b);if(value>0)ratioMin=std::min(ratioMin,value);ratioMax=std::max(ratioMax,value);}pythiaOverFlux->SetMinimum(ratioMin<1.e299?std::max(1.e-6,.5*ratioMin):1.e-3);pythiaOverFlux->SetMaximum(ratioMax>0?std::max(2.,2.*ratioMax):2.);pythiaOverFlux->Draw("HIST");TLine unity(commonLow,1,commonHigh,1);unity.SetLineStyle(2);unity.SetLineColor(kGray+2);unity.DrawClone();heading("Normalized shape ratio");text(.17,.73,"Both spectra normalized on common support",.028);
 save(fluxComparison,"qa_noon_flux_vs_pythia_egamma");

 out.cd();for(auto& h:owned)h->Write();
 TNamed("plotNotes","Original channel-mapping artifacts retained (notably 49/98); no extra flux or detector correction inside this plotting macro; classes select overlapping photon spectra. Closure ratios use correlated entries. Conditional map errors are not plotted; empty energy slices are zero. n00n/PYTHIA comparison projects h2fluxNoon_wPnohad over stored regular b bins and normalizes both shapes on common support.").Write();
 TNamed("fluxComparisonInput",Form("file=%s; histogram=h2fluxNoon_wPnohad; n00n x converted MeV->GeV",fluxInput)).Write();
 for(auto name:{"physicsScope","energyConvention","weightConvention","noonFluxWeightApplied","noonTrackSelection","classDefinition"})if(auto obj=file->Get(name))obj->Write(name);
 out.Close();std::cout<<"Saved "<<figures<<" PDF/PNG figures in "<<outputPath<<" and "<<rootPath<<std::endl;
}
