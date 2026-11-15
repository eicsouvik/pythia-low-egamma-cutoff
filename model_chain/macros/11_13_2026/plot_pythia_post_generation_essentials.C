// Plot the essential post-PYTHIA diagnostics needed for the PYTHIA--nOOn
// interface and for multiplicity-selected correlation studies.
//
// Run from Code/upc2pc/model with:
//   root -l -b -q plot_pythia_post_generation_essentials.C
//
// A different input file or output directory can be supplied with:
//   root -l -b -q \
//     'plot_pythia_post_generation_essentials.C("file.root","output_dir")'

#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TNamed.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

TString MacroDirectory()
{
  TString macroPath = gSystem->UnixPathName(__FILE__);
  TString directory = gSystem->DirName(macroPath);
  if (directory.IsNull()) directory = ".";
  return directory;
}

TString ResolveFromMacro(const TString &path)
{
  if (gSystem->IsAbsoluteFileName(path.Data())) return path;
  return MacroDirectory() + "/" + path;
}

TString ExpandDatePattern(const TString &path)
{
  TString expanded(path);
  if (!expanded.Contains("mm_dd_yy")) return expanded;

  TDatime now;
  expanded.ReplaceAll(
      "mm_dd_yy",
      Form("%02d_%02d_%02d",
           now.GetMonth(), now.GetDay(), now.GetYear() % 100));
  return expanded;
}

template <class TObjectType>
TObjectType *GetRequired(TFile *input, const char *name)
{
  TObjectType *object = dynamic_cast<TObjectType *>(input->Get(name));
  if (!object) {
    throw std::runtime_error(
        Form("Required ROOT object '%s' was not found in %s",
             name, input->GetName()));
  }
  return object;
}

template <class TObjectType>
TObjectType *GetOptional(TFile *input, const char *name)
{
  return dynamic_cast<TObjectType *>(input->Get(name));
}

TH1 *Clone1D(TFile *input, const char *name, const char *cloneName)
{
  TH1 *source = GetRequired<TH1>(input, name);
  if (source->GetDimension() != 1) {
    throw std::runtime_error(Form("'%s' is not one-dimensional", name));
  }
  TH1 *clone = dynamic_cast<TH1 *>(source->Clone(cloneName));
  clone->SetDirectory(nullptr);
  return clone;
}

TH2 *Clone2D(TFile *input, const char *name, const char *cloneName)
{
  TH2 *source = GetRequired<TH2>(input, name);
  TH2 *clone = dynamic_cast<TH2 *>(source->Clone(cloneName));
  clone->SetDirectory(nullptr);
  return clone;
}

void SetAxisStyle(TH1 *histogram)
{
  histogram->SetStats(false);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
  histogram->GetXaxis()->SetTitleSize(0.050);
  histogram->GetYaxis()->SetTitleSize(0.050);
  histogram->GetXaxis()->SetLabelSize(0.042);
  histogram->GetYaxis()->SetLabelSize(0.042);
  histogram->GetXaxis()->SetTitleOffset(1.08);
  histogram->GetYaxis()->SetTitleOffset(1.28);
}

void StyleHistogram(TH1 *histogram, int color, int marker, int lineStyle = 1)
{
  SetAxisStyle(histogram);
  histogram->SetLineColor(color);
  histogram->SetMarkerColor(color);
  histogram->SetMarkerStyle(marker);
  histogram->SetMarkerSize(0.85);
  histogram->SetLineStyle(lineStyle);
  histogram->SetLineWidth(3);
}

void SetPadStyle(bool logY = false, bool logX = false)
{
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.08);
  gPad->SetTicks(1, 1);
  gPad->SetLogx(logX);
  gPad->SetLogy(logY);
}

void SetColzPadStyle(bool logX = false, bool logZ = false)
{
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.09);
  gPad->SetTicks(1, 1);
  gPad->SetLogx(logX);
  gPad->SetLogz(logZ);
}

void AddAnalysisLabel(double x = 0.18, double y = 0.88,
                      const char *extra = nullptr)
{
  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.040);
  label.DrawLatex(x, y, "PYTHIA 8.3 #gamma+p");
  label.SetTextSize(0.033);
  label.DrawLatex(x, y - 0.055, "10#times100 GeV e+p, Q^{2}<0.01 GeV^{2}");
  if (extra && TString(extra).Length() > 0) {
    label.DrawLatex(x, y - 0.105, extra);
  }
}

void NormalizeDensity(TH1 *histogram)
{
  const double integral = histogram->Integral(1, histogram->GetNbinsX(), "width");
  if (integral <= 0.0) {
    std::cerr << "Warning: histogram '" << histogram->GetName()
              << "' has a non-positive width integral; leaving it unscaled."
              << std::endl;
    return;
  }
  histogram->Scale(1.0 / integral);
}

void NormalizeProbability(TH1 *histogram)
{
  const double integral = histogram->Integral(1, histogram->GetNbinsX());
  if (integral <= 0.0) {
    throw std::runtime_error(
        Form("Histogram '%s' has a non-positive integral",
             histogram->GetName()));
  }
  histogram->Scale(1.0 / integral);
}

void NormalizePerEventDensity(TH1 *histogram, double eventCount)
{
  if (eventCount <= 0.0) {
    std::cerr << "Warning: cannot normalize '" << histogram->GetName()
              << "' per event because its class event count is non-positive; "
              << "leaving it unscaled." << std::endl;
    return;
  }
  histogram->Scale(1.0 / eventCount, "width");
}

double PositiveMinimum(const std::vector<TH1 *> &histograms)
{
  double minimum = 1.0e300;
  for (const TH1 *histogram : histograms) {
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
      const double value = histogram->GetBinContent(bin);
      if (value > 0.0) minimum = std::min(minimum, value);
    }
  }
  return minimum < 1.0e299 ? minimum : 1.0e-12;
}

double Maximum(const std::vector<TH1 *> &histograms)
{
  double maximum = 0.0;
  for (const TH1 *histogram : histograms) {
    maximum = std::max(maximum, histogram->GetMaximum());
  }
  return maximum;
}

void SetOverlayRange(TH1 *frame, const std::vector<TH1 *> &histograms,
                     bool logarithmic)
{
  const double maximum = std::max(1.0e-12, Maximum(histograms));
  if (logarithmic) {
    frame->SetMinimum(std::max(1.0e-12, 0.40 * PositiveMinimum(histograms)));
    frame->SetMaximum(8.0 * maximum);
  } else {
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.25 * maximum);
  }
}

TLegend *MakeLegend(double x1, double y1, double x2, double y2)
{
  TLegend *legend = new TLegend(x1, y1, x2, y2);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  return legend;
}

void SaveCanvas(TCanvas *canvas, const TString &outputDirectory,
                const char *stem)
{
  canvas->SaveAs(outputDirectory + "/" + stem + ".pdf");
  canvas->SaveAs(outputDirectory + "/" + stem + ".png");
}

double EventCountInNTofRange(const TH1 *histogram,
                             double minimum, double maximum)
{
  const int firstBin = histogram->GetXaxis()->FindFixBin(minimum + 1.0e-9);
  const int lastBin = histogram->GetXaxis()->FindFixBin(maximum - 1.0e-9);
  return histogram->Integral(firstBin, lastBin);
}

TH1D *ProjectLongRangeDeltaPhi(TH2 *histogram, const char *name)
{
  const double lower = std::max(1.0, histogram->GetXaxis()->GetXmin());
  const double upper = std::min(3.0, histogram->GetXaxis()->GetXmax());
  if (upper <= lower) {
    throw std::runtime_error(
        Form("'%s' does not cover 1<|Delta-eta|<3", histogram->GetName()));
  }
  const int firstBin = histogram->GetXaxis()->FindFixBin(lower + 1.0e-9);
  const int lastBin = histogram->GetXaxis()->FindFixBin(upper - 1.0e-9);
  TH1D *projection = histogram->ProjectionY(name, firstBin, lastBin, "e");
  projection->SetDirectory(nullptr);
  NormalizeDensity(projection);
  return projection;
}

} // namespace

void plot_pythia_post_generation_essentials(
    const char *inputFile =
        "PYTHIA_root_files/"
        "egamma_less1p1gev_ep10x100_test17-18_updated_analysis_macro.root",
    const char *outputDirectory = "PYTHIA_post_generation_plots_mm_dd_yy")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(100);

  const TString inputPath = ResolveFromMacro(inputFile);
  const TString outputPath = ResolveFromMacro(
      ExpandDatePattern(outputDirectory));
  gSystem->mkdir(outputPath, true);

  TFile *input = TFile::Open(inputPath, "READ");
  if (!input || input->IsZombie()) {
    throw std::runtime_error(Form("Could not open %s", inputPath.Data()));
  }

  std::vector<TCanvas *> canvases;

  try {
    // ------------------------------------------------------------------
    // Read all post-PYTHIA response products.
    // ------------------------------------------------------------------
    TH1 *egammaGeneratorAll = Clone1D(
        input, "h_Egamma_generator_all", "plot_Egamma_generator_all");
    TH1 *egammaGeneratorClass0 = Clone1D(
        input, "h_Egamma_generator_class0", "plot_Egamma_generator_class0");
    TH1 *egammaGeneratorClass1 = Clone1D(
        input, "h_Egamma_generator_class1", "plot_Egamma_generator_class1");
    TH1 *egammaGeneratorClass2 = Clone1D(
        input, "h_Egamma_generator_class2", "plot_Egamma_generator_class2");

    TH1 *egammaUpcAll = Clone1D(
        input, "h_Egamma_upc_all", "plot_Egamma_upc_all");
    TH1 *egammaUpcClass0 = Clone1D(
        input, "h_Egamma_upc_class0", "plot_Egamma_upc_class0");
    TH1 *egammaUpcClass1 = Clone1D(
        input, "h_Egamma_upc_class1", "plot_Egamma_upc_class1");
    TH1 *egammaUpcClass2 = Clone1D(
        input, "h_Egamma_upc_class2", "plot_Egamma_upc_class2");

    TH2 *egammaNtofGenerator = Clone2D(
        input, "h2_Egamma_Ntof_generator", "plot_Egamma_Ntof_generator");
    TH2 *egammaNtofUpc = Clone2D(
        input, "h2_Egamma_Ntof_upc", "plot_Egamma_Ntof_upc");

    TH1 *nTofGenerator = Clone1D(
        input, "h_Ntof_generator", "plot_Ntof_generator");
    TH1 *meanNTof = Clone1D(
        input, "p_meanNtof_vs_Egamma", "plot_meanNtof_vs_Egamma");
    TH1 *classProbability0 = Clone1D(
        input, "h_class0_probability_vs_Egamma", "plot_class0_probability");
    TH1 *classProbability1 = Clone1D(
        input, "h_class1_probability_vs_Egamma", "plot_class1_probability");
    TH1 *classProbability2 = Clone1D(
        input, "h_class2_probability_vs_Egamma", "plot_class2_probability");

    TH1 *wGammaPAll = Clone1D(
        input, "h_WgammaP_all", "plot_WgammaP_all");
    TH1 *wGammaPClass0 = Clone1D(
        input, "h_WgammaP_class0", "plot_WgammaP_class0");
    TH1 *wGammaPClass1 = Clone1D(
        input, "h_WgammaP_class1", "plot_WgammaP_class1");

    TH1 *q2 = Clone1D(input, "q2hist", "plot_q2");
    TH1 *pt = Clone1D(input, "h_pt_all", "plot_pt_all");
    TH1 *npartRaw = Clone1D(input, "h_npart_raw", "plot_npart_raw");
    TH1 *npartTof = Clone1D(
        input, "h_npart_tof_matched", "plot_npart_tof_matched");
    TH1 *etaInclusive = Clone1D(input, "h_dNdEta2", "plot_eta_inclusive");
    TH1 *etaClass0 = Clone1D(input, "h_dNdEta2_cent0", "plot_eta_class0");
    TH1 *etaClass1 = Clone1D(input, "h_dNdEta2_cent1", "plot_eta_class1");
    TH1 *etaClass2 = Clone1D(input, "h_dNdEta2_cent2", "plot_eta_class2");

    TH2 *detaDphiClass0 = Clone2D(
        input, "h_deta_dphi_cent0", "plot_deta_dphi_class0");
    TH2 *detaDphiClass1 = Clone2D(
        input, "h_deta_dphi_cent1", "plot_deta_dphi_class1");
    TH2 *detaDphiClass2 = Clone2D(
        input, "h_deta_dphi_cent2", "plot_deta_dphi_class2");

    TGraph *sourceWeights = GetOptional<TGraph>(input, "gWeights");
    TGraph *weightGraph = sourceWeights
        ? dynamic_cast<TGraph *>(sourceWeights->Clone("plot_upc_weights"))
        : nullptr;
    TNamed *weightStatusObject = GetOptional<TNamed>(
        input, "noonFluxWeightApplied");
    const TString weightStatus = weightStatusObject
        ? weightStatusObject->GetTitle() : "status unavailable";

    const int colors[4] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2};
    const int markers[4] = {20, 21, 20, 33};
    const char *labels[4] = {
        "Inclusive", "N_{TOF} #leq 5",
        "5 < N_{TOF} #leq 15", "N_{TOF} > 15"};

    // ------------------------------------------------------------------
    // 1. Photon-flux weighting validation.
    // ------------------------------------------------------------------
    TCanvas *weightCanvas = new TCanvas(
        "c_pythia_weights", "Photon-energy weighting", 1500, 650);
    weightCanvas->Divide(2, 1);

    weightCanvas->cd(1);
    SetPadStyle(false, true);
    if (weightGraph) {
      weightGraph->SetLineColor(kBlue + 1);
      weightGraph->SetLineWidth(3);
      weightGraph->SetMarkerColor(kBlue + 1);
      weightGraph->SetMarkerStyle(20);
      weightGraph->SetTitle(
          ";E_{#gamma} (MeV);nOOn photon-flux weight");
      weightGraph->Draw("AL");
      weightGraph->GetXaxis()->SetLimits(10.0, 5000.0);
    } else {
      TH1F *frame = gPad->DrawFrame(10.0, 0.0, 5000.0, 1.0);
      frame->SetTitle(";E_{#gamma} (MeV);nOOn photon-flux weight");
      TLatex missing;
      missing.SetNDC();
      missing.SetTextAlign(22);
      missing.SetTextSize(0.045);
      missing.DrawLatex(0.52, 0.52, "gWeights is not stored in this file");
    }
    AddAnalysisLabel(0.18, 0.88, Form("Flux weight: %s", weightStatus.Data()));

    weightCanvas->cd(2);
    SetPadStyle(true, true);
    TH1 *egammaGeneratorShape = dynamic_cast<TH1 *>(
        egammaGeneratorAll->Clone("plot_Egamma_generator_shape"));
    TH1 *egammaUpcShape = dynamic_cast<TH1 *>(
        egammaUpcAll->Clone("plot_Egamma_upc_shape"));
    egammaGeneratorShape->SetDirectory(nullptr);
    egammaUpcShape->SetDirectory(nullptr);
    NormalizeDensity(egammaGeneratorShape);
    NormalizeDensity(egammaUpcShape);
    StyleHistogram(egammaGeneratorShape, kBlue + 1, 20);
    StyleHistogram(egammaUpcShape, kRed + 1, 24, 2);
    egammaGeneratorShape->GetXaxis()->SetRangeUser(10.0, 5000.0);
    egammaGeneratorShape->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    egammaGeneratorShape->GetYaxis()->SetTitle("unit-normalized density");
    SetOverlayRange(egammaGeneratorShape,
                    {egammaGeneratorShape, egammaUpcShape}, true);
    egammaGeneratorShape->Draw("HIST");
    egammaUpcShape->Draw("HIST SAME");
    TLegend *weightLegend = MakeLegend(0.49, 0.70, 0.91, 0.87);
    weightLegend->AddEntry(egammaGeneratorShape, "Generator spectrum", "l");
    weightLegend->AddEntry(egammaUpcShape, "After nOOn flux weight", "l");
    weightLegend->Draw();
    AddAnalysisLabel();
    SaveCanvas(weightCanvas, outputPath, "01_photon_energy_weighting");
    canvases.push_back(weightCanvas);

    // ------------------------------------------------------------------
    // 2. E_gamma--N_TOF response matrices before and after UPC weighting.
    // ------------------------------------------------------------------
    TCanvas *responseMatrixCanvas = new TCanvas(
        "c_pythia_response_matrices", "Egamma versus NTOF", 1500, 650);
    responseMatrixCanvas->Divide(2, 1);

    responseMatrixCanvas->cd(1);
    SetColzPadStyle(true, true);
    egammaNtofGenerator->GetXaxis()->SetRangeUser(10.0, 5000.0);
    egammaNtofGenerator->GetYaxis()->SetRangeUser(-0.5, 30.5);
    egammaNtofGenerator->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    egammaNtofGenerator->GetYaxis()->SetTitle("N_{TOF}");
    egammaNtofGenerator->GetZaxis()->SetTitle("generator-weighted events");
    egammaNtofGenerator->Draw("COLZ");
    AddAnalysisLabel(0.17, 0.87, "Generator response");

    responseMatrixCanvas->cd(2);
    SetColzPadStyle(true, true);
    egammaNtofUpc->GetXaxis()->SetRangeUser(10.0, 5000.0);
    egammaNtofUpc->GetYaxis()->SetRangeUser(-0.5, 30.5);
    egammaNtofUpc->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    egammaNtofUpc->GetYaxis()->SetTitle("N_{TOF}");
    egammaNtofUpc->GetZaxis()->SetTitle("UPC-weighted events");
    egammaNtofUpc->Draw("COLZ");
    AddAnalysisLabel(0.17, 0.87, "After nOOn flux weight");
    SaveCanvas(responseMatrixCanvas, outputPath, "02_egamma_ntof_response");
    canvases.push_back(responseMatrixCanvas);

    // ------------------------------------------------------------------
    // 3. Multiplicity response: P(N_TOF), mean N_TOF, and class response.
    // ------------------------------------------------------------------
    TCanvas *multiplicityCanvas = new TCanvas(
        "c_pythia_multiplicity_response", "Multiplicity response", 1800, 600);
    multiplicityCanvas->Divide(3, 1);

    multiplicityCanvas->cd(1);
    SetPadStyle(true, false);
    NormalizeProbability(nTofGenerator);
    StyleHistogram(nTofGenerator, kBlue + 1, 20);
    nTofGenerator->GetXaxis()->SetRangeUser(-0.5, 30.5);
    nTofGenerator->GetXaxis()->SetTitle("N_{TOF}");
    nTofGenerator->GetYaxis()->SetTitle("event probability");
    nTofGenerator->SetMinimum(
        std::max(1.0e-12, 0.4 * PositiveMinimum({nTofGenerator})));
    nTofGenerator->SetMaximum(5.0 * nTofGenerator->GetMaximum());
    nTofGenerator->Draw("E1");
    AddAnalysisLabel(0.19, 0.88, "#epsilon_{TOF}=65% binomial response");

    multiplicityCanvas->cd(2);
    SetPadStyle(false, true);
    StyleHistogram(meanNTof, kBlack, 20);
    meanNTof->GetXaxis()->SetRangeUser(10.0, 5000.0);
    meanNTof->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    meanNTof->GetYaxis()->SetTitle("#LTN_{TOF}#GT");
    meanNTof->SetMinimum(0.0);
    meanNTof->SetMaximum(1.25 * meanNTof->GetMaximum());
    meanNTof->Draw("E1");
    AddAnalysisLabel();

    multiplicityCanvas->cd(3);
    SetPadStyle(false, true);
    StyleHistogram(classProbability0, colors[1], markers[1]);
    StyleHistogram(classProbability1, colors[2], markers[2]);
    StyleHistogram(classProbability2, colors[3], markers[3]);
    classProbability0->GetXaxis()->SetRangeUser(10.0, 5000.0);
    classProbability0->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    classProbability0->GetYaxis()->SetTitle("P(class|E_{#gamma})");
    classProbability0->SetMinimum(0.0);
    classProbability0->SetMaximum(1.08);
    classProbability0->Draw("E1");
    classProbability1->Draw("E1 SAME");
    classProbability2->Draw("E1 SAME");
    TLegend *probabilityLegend = MakeLegend(0.43, 0.65, 0.92, 0.88);
    probabilityLegend->AddEntry(classProbability0, labels[1], "lep");
    probabilityLegend->AddEntry(classProbability1, labels[2], "lep");
    probabilityLegend->AddEntry(classProbability2, labels[3], "lep");
    probabilityLegend->Draw();
    SaveCanvas(multiplicityCanvas, outputPath, "03_multiplicity_response");
    canvases.push_back(multiplicityCanvas);

    // ------------------------------------------------------------------
    // 4. Conditional E_gamma spectra: rates and normalized shapes.
    // ------------------------------------------------------------------
    TCanvas *egammaClassCanvas = new TCanvas(
        "c_pythia_egamma_classes", "Egamma by multiplicity", 1500, 1200);
    egammaClassCanvas->Divide(2, 2);

    std::vector<TH1 *> generatorSpectra = {
        egammaGeneratorAll, egammaGeneratorClass0,
        egammaGeneratorClass1, egammaGeneratorClass2};
    std::vector<TH1 *> upcSpectra = {
        egammaUpcAll, egammaUpcClass0, egammaUpcClass1, egammaUpcClass2};
    for (int index = 0; index < 4; ++index) {
      StyleHistogram(generatorSpectra[index], colors[index], markers[index]);
      StyleHistogram(upcSpectra[index], colors[index], markers[index]);
    }

    egammaClassCanvas->cd(1);
    SetPadStyle(true, true);
    egammaGeneratorAll->GetXaxis()->SetRangeUser(10.0, 5000.0);
    egammaGeneratorAll->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    egammaGeneratorAll->GetYaxis()->SetTitle("generator-weighted events");
    SetOverlayRange(egammaGeneratorAll, generatorSpectra, true);
    egammaGeneratorAll->Draw("HIST");
    for (int index = 1; index < 4; ++index)
      generatorSpectra[index]->Draw("HIST SAME");
    TLegend *generatorLegend = MakeLegend(0.48, 0.61, 0.91, 0.88);
    for (int index = 0; index < 4; ++index)
      generatorLegend->AddEntry(generatorSpectra[index], labels[index], "l");
    generatorLegend->Draw();
    AddAnalysisLabel(0.18, 0.88, "Before UPC weighting");

    egammaClassCanvas->cd(2);
    SetPadStyle(true, true);
    egammaUpcAll->GetXaxis()->SetRangeUser(10.0, 5000.0);
    egammaUpcAll->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    egammaUpcAll->GetYaxis()->SetTitle("UPC-weighted events");
    SetOverlayRange(egammaUpcAll, upcSpectra, true);
    egammaUpcAll->Draw("HIST");
    for (int index = 1; index < 4; ++index)
      upcSpectra[index]->Draw("HIST SAME");
    TLegend *upcLegend = MakeLegend(0.48, 0.61, 0.91, 0.88);
    for (int index = 0; index < 4; ++index)
      upcLegend->AddEntry(upcSpectra[index], labels[index], "l");
    upcLegend->Draw();
    AddAnalysisLabel(0.18, 0.88, "After nOOn flux weight");

    std::vector<TH1 *> generatorShapes;
    std::vector<TH1 *> upcShapes;
    for (int index = 0; index < 4; ++index) {
      TH1 *generatorShape = dynamic_cast<TH1 *>(generatorSpectra[index]->Clone(
          Form("plot_generator_Egamma_shape_%d", index)));
      TH1 *upcShape = dynamic_cast<TH1 *>(upcSpectra[index]->Clone(
          Form("plot_upc_Egamma_shape_%d", index)));
      generatorShape->SetDirectory(nullptr);
      upcShape->SetDirectory(nullptr);
      NormalizeDensity(generatorShape);
      NormalizeDensity(upcShape);
      StyleHistogram(generatorShape, colors[index], markers[index]);
      StyleHistogram(upcShape, colors[index], markers[index]);
      generatorShapes.push_back(generatorShape);
      upcShapes.push_back(upcShape);
    }

    egammaClassCanvas->cd(3);
    SetPadStyle(true, true);
    generatorShapes[0]->GetXaxis()->SetRangeUser(10.0, 5000.0);
    generatorShapes[0]->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    generatorShapes[0]->GetYaxis()->SetTitle("unit-normalized density");
    SetOverlayRange(generatorShapes[0], generatorShapes, true);
    generatorShapes[0]->Draw("HIST");
    for (int index = 1; index < 4; ++index)
      generatorShapes[index]->Draw("HIST SAME");
    AddAnalysisLabel(0.18, 0.88, "Generator-level class shapes");

    egammaClassCanvas->cd(4);
    SetPadStyle(true, true);
    upcShapes[0]->GetXaxis()->SetRangeUser(10.0, 5000.0);
    upcShapes[0]->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    upcShapes[0]->GetYaxis()->SetTitle("unit-normalized density");
    SetOverlayRange(upcShapes[0], upcShapes, true);
    upcShapes[0]->Draw("HIST");
    for (int index = 1; index < 4; ++index)
      upcShapes[index]->Draw("HIST SAME");
    AddAnalysisLabel(0.18, 0.88, "UPC-weighted class shapes");
    SaveCanvas(egammaClassCanvas, outputPath, "04_egamma_by_multiplicity");
    canvases.push_back(egammaClassCanvas);

    // ------------------------------------------------------------------
    // 5. W_gamma-p distributions for the two analysis multiplicity bins.
    // ------------------------------------------------------------------
    TCanvas *wCanvas = new TCanvas(
        "c_pythia_wgammap", "W gamma-p", 1500, 650);
    wCanvas->Divide(2, 1);
    StyleHistogram(wGammaPAll, colors[0], markers[0]);
    StyleHistogram(wGammaPClass0, colors[1], markers[1]);
    StyleHistogram(wGammaPClass1, colors[2], markers[2]);

    wCanvas->cd(1);
    SetPadStyle(true, false);
    wGammaPAll->GetXaxis()->SetTitle("W_{#gamma p} (GeV)");
    wGammaPAll->GetYaxis()->SetTitle("UPC-weighted events");
    SetOverlayRange(wGammaPAll,
                    {wGammaPAll, wGammaPClass0, wGammaPClass1}, true);
    wGammaPAll->Draw("HIST");
    wGammaPClass0->Draw("HIST SAME");
    wGammaPClass1->Draw("HIST SAME");
    TLegend *wLegend = MakeLegend(0.52, 0.65, 0.91, 0.87);
    wLegend->AddEntry(wGammaPAll, labels[0], "l");
    wLegend->AddEntry(wGammaPClass0, labels[1], "l");
    wLegend->AddEntry(wGammaPClass1, labels[2], "l");
    wLegend->Draw();
    AddAnalysisLabel();

    wCanvas->cd(2);
    SetPadStyle(false, false);
    TH1 *wAllShape = dynamic_cast<TH1 *>(wGammaPAll->Clone("plot_W_all_shape"));
    TH1 *wClass0Shape = dynamic_cast<TH1 *>(
        wGammaPClass0->Clone("plot_W_class0_shape"));
    TH1 *wClass1Shape = dynamic_cast<TH1 *>(
        wGammaPClass1->Clone("plot_W_class1_shape"));
    wAllShape->SetDirectory(nullptr);
    wClass0Shape->SetDirectory(nullptr);
    wClass1Shape->SetDirectory(nullptr);
    NormalizeDensity(wAllShape);
    NormalizeDensity(wClass0Shape);
    NormalizeDensity(wClass1Shape);
    StyleHistogram(wAllShape, colors[0], markers[0]);
    StyleHistogram(wClass0Shape, colors[1], markers[1]);
    StyleHistogram(wClass1Shape, colors[2], markers[2]);
    wAllShape->GetXaxis()->SetTitle("W_{#gamma p} (GeV)");
    wAllShape->GetYaxis()->SetTitle("unit-normalized density");
    SetOverlayRange(wAllShape, {wAllShape, wClass0Shape, wClass1Shape}, false);
    wAllShape->Draw("HIST");
    wClass0Shape->Draw("HIST SAME");
    wClass1Shape->Draw("HIST SAME");
    wLegend->Draw();
    AddAnalysisLabel();
    SaveCanvas(wCanvas, outputPath, "05_wgammap_by_multiplicity");
    canvases.push_back(wCanvas);

    // ------------------------------------------------------------------
    // 6. Basic PYTHIA final-state and detector-response QA.
    // ------------------------------------------------------------------
    // nTofGenerator was normalized above. Obtain the original class counts
    // directly from the response matrix Y projection to retain normalization.
    TH1D *nTofCounts = egammaNtofGenerator->ProjectionY(
        "plot_Ntof_counts_from_response", 1, egammaNtofGenerator->GetNbinsX(), "e");
    nTofCounts->SetDirectory(nullptr);
    const double eventCount0 = EventCountInNTofRange(nTofCounts, 0.0, 5.0);
    const double eventCount1 = EventCountInNTofRange(nTofCounts, 6.0, 15.0);
    const double eventCount2 = EventCountInNTofRange(
        nTofCounts, 16.0, nTofCounts->GetXaxis()->GetXmax());
    const double eventCountAll = nTofCounts->Integral(1, nTofCounts->GetNbinsX());

    TCanvas *qaCanvas = new TCanvas(
        "c_pythia_final_state_qa", "PYTHIA final-state QA", 1500, 1200);
    qaCanvas->Divide(2, 2);

    qaCanvas->cd(1);
    SetPadStyle(true, false);
    NormalizeDensity(q2);
    StyleHistogram(q2, kBlue + 1, 20);
    q2->GetXaxis()->SetRangeUser(0.0, 0.01);
    q2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    q2->GetYaxis()->SetTitle("unit-normalized density");
    q2->SetMinimum(std::max(1.0e-12, 0.4 * PositiveMinimum({q2})));
    q2->SetMaximum(5.0 * q2->GetMaximum());
    q2->Draw("HIST");
    AddAnalysisLabel();

    qaCanvas->cd(2);
    SetPadStyle(true, false);
    NormalizePerEventDensity(pt, eventCountAll);
    StyleHistogram(pt, kRed + 1, 20);
    pt->GetXaxis()->SetRangeUser(0.0, 5.0);
    pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    pt->GetYaxis()->SetTitle("1/N_{evt} dN/dp_{T}");
    pt->SetMinimum(std::max(1.0e-12, 0.4 * PositiveMinimum({pt})));
    pt->SetMaximum(5.0 * pt->GetMaximum());
    pt->Draw("HIST");
    AddAnalysisLabel();

    qaCanvas->cd(3);
    SetPadStyle(true, false);
    NormalizeProbability(npartRaw);
    NormalizeProbability(npartTof);
    StyleHistogram(npartRaw, kBlack, 20);
    StyleHistogram(npartTof, kBlue + 1, 24);
    npartRaw->GetXaxis()->SetRangeUser(-0.5, 50.5);
    npartRaw->GetXaxis()->SetTitle("multiplicity");
    npartRaw->GetYaxis()->SetTitle("event probability");
    SetOverlayRange(npartRaw, {npartRaw, npartTof}, true);
    npartRaw->Draw("HIST");
    npartTof->Draw("HIST SAME");
    TLegend *multiplicityLegend = MakeLegend(0.51, 0.70, 0.91, 0.87);
    multiplicityLegend->AddEntry(npartRaw, "Raw final-state multiplicity", "l");
    multiplicityLegend->AddEntry(npartTof, "TOF-matched multiplicity", "l");
    multiplicityLegend->Draw();
    AddAnalysisLabel();

    qaCanvas->cd(4);
    SetPadStyle(false, false);
    NormalizePerEventDensity(etaInclusive, eventCountAll);
    NormalizePerEventDensity(etaClass0, eventCount0);
    NormalizePerEventDensity(etaClass1, eventCount1);
    NormalizePerEventDensity(etaClass2, eventCount2);
    std::vector<TH1 *> etaHistograms = {
        etaInclusive, etaClass0, etaClass1, etaClass2};
    for (int index = 0; index < 4; ++index)
      StyleHistogram(etaHistograms[index], colors[index], markers[index]);
    etaInclusive->GetXaxis()->SetRangeUser(-1.05, 1.05);
    etaInclusive->GetXaxis()->SetTitle("#eta");
    etaInclusive->GetYaxis()->SetTitle("1/N_{evt} dN/d#eta");
    SetOverlayRange(etaInclusive, etaHistograms, false);
    etaInclusive->Draw("HIST");
    for (int index = 1; index < 4; ++index)
      etaHistograms[index]->Draw("HIST SAME");
    TLegend *etaLegend = MakeLegend(0.51, 0.60, 0.91, 0.87);
    for (int index = 0; index < 4; ++index)
      etaLegend->AddEntry(etaHistograms[index], labels[index], "l");
    etaLegend->Draw();
    AddAnalysisLabel(0.18, 0.88, "|#eta|<1, 0.2<p_{T}<2 GeV/c");
    SaveCanvas(qaCanvas, outputPath, "06_final_state_and_detector_QA");
    canvases.push_back(qaCanvas);

    // ------------------------------------------------------------------
    // 7. Delta-eta--Delta-phi maps for all multiplicity classes.
    // ------------------------------------------------------------------
    TCanvas *correlationCanvas = new TCanvas(
        "c_pythia_correlation_maps", "Deta-Dphi maps", 1800, 600);
    correlationCanvas->Divide(3, 1);
    TH2 *correlationHistograms[3] = {
        detaDphiClass0, detaDphiClass1, detaDphiClass2};
    for (int index = 0; index < 3; ++index) {
      correlationCanvas->cd(index + 1);
      SetColzPadStyle(false, false);
      correlationHistograms[index]->GetXaxis()->SetRangeUser(0.0, 3.0);
      correlationHistograms[index]->GetXaxis()->SetTitle("|#Delta#eta|");
      correlationHistograms[index]->GetYaxis()->SetTitle("#Delta#varphi");
      correlationHistograms[index]->GetZaxis()->SetTitle("weighted pairs");
      correlationHistograms[index]->Draw("COLZ");
      AddAnalysisLabel(0.17, 0.87, labels[index + 1]);
    }
    SaveCanvas(correlationCanvas, outputPath, "07_deta_dphi_by_multiplicity");
    canvases.push_back(correlationCanvas);

    // ------------------------------------------------------------------
    // 8. Long-range Delta-phi projections, 1<|Delta-eta|<3.
    // ------------------------------------------------------------------
    TH1D *deltaPhiClass0 = ProjectLongRangeDeltaPhi(
        detaDphiClass0, "plot_delta_phi_class0");
    TH1D *deltaPhiClass1 = ProjectLongRangeDeltaPhi(
        detaDphiClass1, "plot_delta_phi_class1");
    TH1D *deltaPhiClass2 = ProjectLongRangeDeltaPhi(
        detaDphiClass2, "plot_delta_phi_class2");
    StyleHistogram(deltaPhiClass0, colors[1], markers[1]);
    StyleHistogram(deltaPhiClass1, colors[2], markers[2]);
    StyleHistogram(deltaPhiClass2, colors[3], markers[3]);

    TCanvas *deltaPhiCanvas = new TCanvas(
        "c_pythia_long_range_dphi", "Long-range Delta phi", 850, 700);
    SetPadStyle(false, false);
    deltaPhiClass0->GetXaxis()->SetTitle("#Delta#varphi");
    deltaPhiClass0->GetYaxis()->SetTitle("unit-normalized pair density");
    SetOverlayRange(deltaPhiClass0,
                    {deltaPhiClass0, deltaPhiClass1, deltaPhiClass2}, false);
    deltaPhiClass0->Draw("E1");
    deltaPhiClass1->Draw("E1 SAME");
    deltaPhiClass2->Draw("E1 SAME");
    TLegend *deltaPhiLegend = MakeLegend(0.50, 0.65, 0.91, 0.87);
    deltaPhiLegend->AddEntry(deltaPhiClass0, labels[1], "lep");
    deltaPhiLegend->AddEntry(deltaPhiClass1, labels[2], "lep");
    deltaPhiLegend->AddEntry(deltaPhiClass2, labels[3], "lep");
    deltaPhiLegend->Draw();
    AddAnalysisLabel(0.18, 0.88, "1<|#Delta#eta|<3");
    SaveCanvas(deltaPhiCanvas, outputPath, "08_long_range_dphi");
    canvases.push_back(deltaPhiCanvas);

    // Write a convenient multi-page PDF containing every canvas.
    const TString combinedPdf = outputPath + "/PYTHIA_essential_plots.pdf";
    for (std::size_t index = 0; index < canvases.size(); ++index) {
      TString destination = combinedPdf;
      if (index == 0) destination += "(";
      if (index + 1 == canvases.size()) destination += ")";
      canvases[index]->Print(destination);
    }

    std::cout << "\nCreated PYTHIA post-generation plots in: "
              << outputPath << std::endl;
    std::cout << "Combined PDF: " << combinedPdf << std::endl;
    std::cout << "Input file: " << inputPath << std::endl;
    std::cout << "nOOn flux-weight status: " << weightStatus << std::endl;

  } catch (const std::exception &error) {
    std::cerr << "plot_pythia_post_generation_essentials failed: "
              << error.what() << std::endl;
  }

  input->Close();
}
