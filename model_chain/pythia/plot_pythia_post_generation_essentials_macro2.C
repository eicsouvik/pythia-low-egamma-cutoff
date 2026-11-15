// Extended post-generation QA for the November 14 PYTHIA production.
//
// Run from Code/upc2pc/model with:
//   root -l -b -q plot_pythia_post_generation_essentials_macro2.C
//
// Optional arguments:
//   root -l -b -q \
//     'plot_pythia_post_generation_essentials_macro2.C(
//        "PYTHIA_root_files/input.root",
//        "n00n_flux_from_Ashik_da/noon_flux.root",
//        "output_folder_mm_dd_yy")'

// Reuse the established essential-plot style and the original eight canvases.
#include "plot_pythia_post_generation_essentials.C"

#include "TLine.h"

#include <limits>

namespace {

void ZeroOutside(TH1 *histogram, int firstBin, int lastBin)
{
  for (int bin = 0; bin <= histogram->GetNbinsX() + 1; ++bin) {
    if (bin < firstBin || bin > lastBin) {
      histogram->SetBinContent(bin, 0.0);
      histogram->SetBinError(bin, 0.0);
    }
  }
}

void NormalizeCountsAsDensity(TH1 *histogram, int firstBin, int lastBin)
{
  ZeroOutside(histogram, firstBin, lastBin);
  const double countIntegral = histogram->Integral(firstBin, lastBin);
  if (!(countIntegral > 0.0)) {
    throw std::runtime_error(
        Form("Histogram '%s' has no counts in the common support",
             histogram->GetName()));
  }
  // The stored PYTHIA contents are counts per variable-width bin.
  histogram->Scale(1.0 / countIntegral, "width");
}

void NormalizeStoredDensity(TH1 *histogram, int firstBin, int lastBin)
{
  ZeroOutside(histogram, firstBin, lastBin);
  const double densityIntegral =
      histogram->Integral(firstBin, lastBin, "width");
  if (!(densityIntegral > 0.0)) {
    throw std::runtime_error(
        Form("Density '%s' has no support in the common range",
             histogram->GetName()));
  }
  // The n00n projection is stored as a density in E_gamma.
  histogram->Scale(1.0 / densityIntegral);
}

void AverageDensity(const TH1 *source, double low, double high,
                    double &average, double &error)
{
  double integral = 0.0;
  double error2 = 0.0;
  for (int bin = 1; bin <= source->GetNbinsX(); ++bin) {
    const double sourceLow = source->GetXaxis()->GetBinLowEdge(bin);
    const double sourceHigh = source->GetXaxis()->GetBinUpEdge(bin);
    const double overlap = std::max(
        0.0, std::min(high, sourceHigh) - std::max(low, sourceLow));
    if (overlap <= 0.0) continue;
    integral += source->GetBinContent(bin) * overlap;
    const double term = source->GetBinError(bin) * overlap;
    error2 += term * term;
  }
  const double width = high - low;
  average = width > 0.0 ? integral / width : 0.0;
  error = width > 0.0 ? std::sqrt(error2) / width : 0.0;
}

TH1D *MapDensityToReferenceAxis(const TH1 *sourceDensity,
                                const TH1 *reference,
                                const char *name)
{
  const TAxis *axis = reference->GetXaxis();
  TH1D *mapped = nullptr;
  if (axis->GetXbins()->GetSize() > 0) {
    mapped = new TH1D(
        name, "", axis->GetNbins(), axis->GetXbins()->GetArray());
  } else {
    mapped = new TH1D(
        name, "", axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
  }
  mapped->SetDirectory(nullptr);
  mapped->Sumw2();
  for (int bin = 1; bin <= mapped->GetNbinsX(); ++bin) {
    const double low = mapped->GetXaxis()->GetBinLowEdge(bin);
    const double high = mapped->GetXaxis()->GetBinUpEdge(bin);
    double average = 0.0;
    double error = 0.0;
    AverageDensity(sourceDensity, low, high, average, error);
    mapped->SetBinContent(bin, average);
    mapped->SetBinError(bin, error);
  }
  return mapped;
}

void FindCommonSupport(const TH1 *pythia, const TH1 *target,
                       int &firstBin, int &lastBin,
                       double &minimum, double &maximum)
{
  firstBin = -1;
  lastBin = -1;
  const double targetMinimum = target->GetXaxis()->GetXmin();
  const double targetMaximum = target->GetXaxis()->GetXmax();
  for (int bin = 1; bin <= pythia->GetNbinsX(); ++bin) {
    const double low = pythia->GetXaxis()->GetBinLowEdge(bin);
    const double high = pythia->GetXaxis()->GetBinUpEdge(bin);
    if (pythia->GetBinContent(bin) <= 0.0 ||
        high <= targetMinimum || low >= targetMaximum ||
        low >= 5000.0) continue;
    if (firstBin < 0) firstBin = bin;
    lastBin = bin;
  }
  if (firstBin < 0 || lastBin < firstBin) {
    throw std::runtime_error(
        "PYTHIA and n00n E_gamma spectra have no common positive support");
  }
  minimum = std::max(
      targetMinimum, pythia->GetXaxis()->GetBinLowEdge(firstBin));
  maximum = std::min(
      5000.0,
      std::min(targetMaximum,
               pythia->GetXaxis()->GetBinUpEdge(lastBin)));
}

TH1D *MakeRatio(const TH1 *numerator, const TH1 *denominator,
                const char *name)
{
  TH1D *ratio = dynamic_cast<TH1D *>(numerator->Clone(name));
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICES");
  for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
    const double denominatorValue = denominator->GetBinContent(bin);
    if (denominatorValue <= 0.0) continue;
    const double numeratorValue = numerator->GetBinContent(bin);
    ratio->SetBinContent(bin, numeratorValue / denominatorValue);
    const double numeratorError = numerator->GetBinError(bin);
    const double denominatorError = denominator->GetBinError(bin);
    double relative2 = 0.0;
    if (numeratorValue > 0.0)
      relative2 += std::pow(numeratorError / numeratorValue, 2);
    relative2 += std::pow(denominatorError / denominatorValue, 2);
    ratio->SetBinError(
        bin, numeratorValue / denominatorValue * std::sqrt(relative2));
  }
  return ratio;
}

void DrawStageMultiplicityPanel(TH1 *raw, TH1 *accepted, TH1 *tpc, TH1 *tof,
                                const char *extra)
{
  NormalizeProbability(raw);
  NormalizeProbability(accepted);
  NormalizeProbability(tpc);
  NormalizeProbability(tof);
  StyleHistogram(raw, kBlack, 20);
  StyleHistogram(accepted, kBlue + 1, 21);
  StyleHistogram(tpc, kGreen + 2, 33);
  StyleHistogram(tof, kRed + 1, 24);
  raw->GetXaxis()->SetRangeUser(-0.5, 50.5);
  raw->GetXaxis()->SetTitle("multiplicity");
  raw->GetYaxis()->SetTitle("event probability");
  SetOverlayRange(raw, {raw, accepted, tpc, tof}, true);
  raw->Draw("HIST");
  accepted->Draw("HIST SAME");
  tpc->Draw("HIST SAME");
  tof->Draw("HIST SAME");
  TLegend *legend = MakeLegend(0.51, 0.64, 0.92, 0.88);
  legend->AddEntry(raw, "Raw charged", "l");
  legend->AddEntry(accepted, "Kinematic acceptance", "l");
  legend->AddEntry(tpc, "TPC reconstructed (95%)", "l");
  legend->AddEntry(tof, "TOF matched (65% of TPC)", "l");
  legend->Draw();
  AddAnalysisLabel(0.18, 0.88, extra);
}

void DrawStageEgammaPanel(std::vector<TH1 *> &histograms,
                          const char *extra)
{
  const int colors[4] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
  const int markers[4] = {20, 21, 33, 24};
  const char *labels[4] = {
      "Raw charged > 0", "Accepted > 0",
      "TPC reconstructed > 0", "TOF matched > 0"};
  for (int index = 0; index < 4; ++index) {
    const double integral = histograms[index]->Integral();
    if (integral > 0.0)
      histograms[index]->Scale(1.0 / integral, "width");
    StyleHistogram(histograms[index], colors[index], markers[index]);
    histograms[index]->GetXaxis()->SetRangeUser(10.0, 5000.0);
  }
  histograms[0]->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  histograms[0]->GetYaxis()->SetTitle("unit-normalized density (MeV^{-1})");
  SetOverlayRange(histograms[0], histograms, true);
  histograms[0]->Draw("HIST");
  for (int index = 1; index < 4; ++index)
    histograms[index]->Draw("HIST SAME");
  TLegend *legend = MakeLegend(0.48, 0.61, 0.92, 0.88);
  for (int index = 0; index < 4; ++index)
    legend->AddEntry(histograms[index], labels[index], "l");
  legend->Draw();
  AddAnalysisLabel(0.18, 0.88, extra);
}

void DrawResponseMatrix(TH2 *histogram, const char *extra,
                        const char *zTitle)
{
  SetColzPadStyle(true, true);
  histogram->GetXaxis()->SetRangeUser(10.0, 5000.0);
  histogram->GetYaxis()->SetRangeUser(-0.5, 60.5);
  histogram->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  histogram->GetZaxis()->SetTitle(zTitle);
  histogram->Draw("COLZ");
  AddAnalysisLabel(0.17, 0.87, extra);
}

} // namespace

void plot_pythia_post_generation_essentials_macro2(
    const char *inputFile =
        "PYTHIA_root_files/"
        "egamma_less5gev_ep10x100_test17-18_updated_analysis_macro2.root",
    const char *noonFluxFile = "n00n_flux_from_Ashik_da/noon_flux.root",
    const char *outputDirectory =
        "PYTHIA_post_generation_plots_macro2_mm_dd_yy")
{
  // First produce the established eight essential canvases.
  plot_pythia_post_generation_essentials(inputFile, outputDirectory);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(100);

  const TString inputPath = ResolveFromMacro(inputFile);
  const TString fluxPath = ResolveFromMacro(noonFluxFile);
  const TString outputPath = ResolveFromMacro(
      ExpandDatePattern(outputDirectory));
  gSystem->mkdir(outputPath, true);

  TFile *input = TFile::Open(inputPath, "READ");
  if (!input || input->IsZombie()) {
    throw std::runtime_error(Form("Could not open %s", inputPath.Data()));
  }
  TFile *fluxInput = TFile::Open(fluxPath, "READ");
  if (!fluxInput || fluxInput->IsZombie()) {
    input->Close();
    throw std::runtime_error(Form("Could not open %s", fluxPath.Data()));
  }

  std::vector<TCanvas *> extraCanvases;
  std::vector<TH1 *> outputHistograms;
  std::vector<TH2 *> outputMatrices;

  try {
    // ------------------------------------------------------------------
    // Read multiplicity distributions at every detector-response stage.
    // ------------------------------------------------------------------
    TH1 *nRawGenerator = Clone1D(
        input, "h_Nraw_generator", "plot2_Nraw_generator");
    TH1 *nRawUpc = Clone1D(
        input, "h_Nraw_upc", "plot2_Nraw_upc");
    TH1 *nAcceptGenerator = Clone1D(
        input, "h_Naccept_generator", "plot2_Naccept_generator");
    TH1 *nAcceptUpc = Clone1D(
        input, "h_Naccept_upc", "plot2_Naccept_upc");
    TH1 *nTpcGenerator = Clone1D(
        input, "h_Ntpc_generator", "plot2_Ntpc_generator");
    TH1 *nTpcUpc = Clone1D(
        input, "h_Ntpc_upc", "plot2_Ntpc_upc");
    TH1 *nTofGenerator = Clone1D(
        input, "h_Ntof_generator", "plot2_Ntof_generator");
    TH1 *nTofUpc = Clone1D(
        input, "h_Ntof_upc", "plot2_Ntof_upc");

    outputHistograms.insert(outputHistograms.end(), {
        nRawGenerator, nRawUpc, nAcceptGenerator, nAcceptUpc,
        nTpcGenerator, nTpcUpc, nTofGenerator, nTofUpc});

    TCanvas *stageMultiplicityCanvas = new TCanvas(
        "c_detector_stage_multiplicity",
        "Multiplicity through detector response", 1500, 650);
    stageMultiplicityCanvas->Divide(2, 1);
    stageMultiplicityCanvas->cd(1);
    SetPadStyle(true, false);
    DrawStageMultiplicityPanel(
        nRawGenerator, nAcceptGenerator, nTpcGenerator, nTofGenerator,
        "Generator-weighted response");
    stageMultiplicityCanvas->cd(2);
    SetPadStyle(true, false);
    DrawStageMultiplicityPanel(
        nRawUpc, nAcceptUpc, nTpcUpc, nTofUpc,
        "After n00n flux weight");
    SaveCanvas(stageMultiplicityCanvas, outputPath,
               "09_detector_stage_multiplicity");
    extraCanvases.push_back(stageMultiplicityCanvas);

    // ------------------------------------------------------------------
    // E_gamma distributions after requiring >=1 object at each stage.
    // ------------------------------------------------------------------
    const char *stageNames[4] = {"raw", "accept", "tpc", "tof"};
    std::vector<TH1 *> egammaStageGenerator;
    std::vector<TH1 *> egammaStageUpc;
    for (int stage = 0; stage < 4; ++stage) {
      TH1 *generator = Clone1D(
          input, Form("h_Egamma_%s_generator", stageNames[stage]),
          Form("plot2_Egamma_%s_generator", stageNames[stage]));
      TH1 *upc = Clone1D(
          input, Form("h_Egamma_%s_upc", stageNames[stage]),
          Form("plot2_Egamma_%s_upc", stageNames[stage]));
      egammaStageGenerator.push_back(generator);
      egammaStageUpc.push_back(upc);
      outputHistograms.push_back(generator);
      outputHistograms.push_back(upc);
    }

    TCanvas *stageEgammaCanvas = new TCanvas(
        "c_detector_stage_egamma",
        "Egamma through detector response", 1500, 650);
    stageEgammaCanvas->Divide(2, 1);
    stageEgammaCanvas->cd(1);
    SetPadStyle(true, true);
    DrawStageEgammaPanel(
        egammaStageGenerator, "Generator-weighted response");
    stageEgammaCanvas->cd(2);
    SetPadStyle(true, true);
    DrawStageEgammaPanel(egammaStageUpc, "After n00n flux weight");
    SaveCanvas(stageEgammaCanvas, outputPath,
               "10_detector_stage_egamma");
    extraCanvases.push_back(stageEgammaCanvas);

    // ------------------------------------------------------------------
    // E_gamma-versus-multiplicity matrices at all four response stages.
    // ------------------------------------------------------------------
    const char *matrixStageName[4] = {"Nraw", "Naccept", "Ntpc", "Ntof"};
    const char *matrixStageLabel[4] = {
        "Raw charged", "Kinematic acceptance",
        "TPC reconstructed (95%)", "TOF matched (65% of TPC)"};
    TH2 *generatorMatrices[4];
    TH2 *upcMatrices[4];
    for (int stage = 0; stage < 4; ++stage) {
      generatorMatrices[stage] = Clone2D(
          input, Form("h2_Egamma_%s_generator", matrixStageName[stage]),
          Form("plot2_Egamma_%s_generator", matrixStageName[stage]));
      upcMatrices[stage] = Clone2D(
          input, Form("h2_Egamma_%s_upc", matrixStageName[stage]),
          Form("plot2_Egamma_%s_upc", matrixStageName[stage]));
      outputMatrices.push_back(generatorMatrices[stage]);
      outputMatrices.push_back(upcMatrices[stage]);
    }

    TCanvas *generatorMatrixCanvas = new TCanvas(
        "c_detector_matrices_generator",
        "Generator detector-response matrices", 1600, 1200);
    generatorMatrixCanvas->Divide(2, 2);
    for (int stage = 0; stage < 4; ++stage) {
      generatorMatrixCanvas->cd(stage + 1);
      DrawResponseMatrix(
          generatorMatrices[stage], matrixStageLabel[stage],
          "generator-weighted events");
    }
    SaveCanvas(generatorMatrixCanvas, outputPath,
               "11_egamma_multiplicity_matrices_generator");
    extraCanvases.push_back(generatorMatrixCanvas);

    TCanvas *upcMatrixCanvas = new TCanvas(
        "c_detector_matrices_upc",
        "UPC detector-response matrices", 1600, 1200);
    upcMatrixCanvas->Divide(2, 2);
    for (int stage = 0; stage < 4; ++stage) {
      upcMatrixCanvas->cd(stage + 1);
      DrawResponseMatrix(
          upcMatrices[stage], matrixStageLabel[stage],
          "n00n-weighted events");
    }
    SaveCanvas(upcMatrixCanvas, outputPath,
               "12_egamma_multiplicity_matrices_upc");
    extraCanvases.push_back(upcMatrixCanvas);

    // ------------------------------------------------------------------
    // Independent n00n photon-flux projection versus PYTHIA input E_gamma.
    // ------------------------------------------------------------------
    TH1 *pythiaInput = Clone1D(
        input, "h_Egamma_generator_all", "plot2_pythia_input_density");
    TH2 *noonFlux2D = GetRequired<TH2>(
        fluxInput, "h2fluxNoon_wPnohad");
    TH1D *noonProjection = noonFlux2D->ProjectionX(
        "plot2_noon_flux_projection", 1, noonFlux2D->GetNbinsY(), "e");
    noonProjection->SetDirectory(nullptr);

    int firstCommonBin = -1;
    int lastCommonBin = -1;
    double commonMinimum = 0.0;
    double commonMaximum = 0.0;
    FindCommonSupport(
        pythiaInput, noonProjection, firstCommonBin, lastCommonBin,
        commonMinimum, commonMaximum);
    TH1D *noonOnPythiaAxis = MapDensityToReferenceAxis(
        noonProjection, pythiaInput, "plot2_noon_density_on_pythia_axis");
    NormalizeCountsAsDensity(
        pythiaInput, firstCommonBin, lastCommonBin);
    NormalizeStoredDensity(
        noonOnPythiaAxis, firstCommonBin, lastCommonBin);
    TH1D *pythiaOverNoon = MakeRatio(
        pythiaInput, noonOnPythiaAxis, "plot2_pythia_over_noon");

    outputHistograms.push_back(pythiaInput);
    outputHistograms.push_back(noonProjection);
    outputHistograms.push_back(noonOnPythiaAxis);
    outputHistograms.push_back(pythiaOverNoon);

    TCanvas *fluxComparisonCanvas = new TCanvas(
        "c_noon_pythia_egamma_comparison",
        "n00n versus PYTHIA photon energy", 1500, 650);
    fluxComparisonCanvas->Divide(2, 1);
    fluxComparisonCanvas->cd(1);
    SetPadStyle(true, true);
    StyleHistogram(pythiaInput, kBlue + 1, 20);
    StyleHistogram(noonOnPythiaAxis, kBlack, 24, 2);
    pythiaInput->GetXaxis()->SetRangeUser(commonMinimum, commonMaximum);
    pythiaInput->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    pythiaInput->GetYaxis()->SetTitle(
        "unit-normalized density (MeV^{-1})");
    SetOverlayRange(
        pythiaInput, {pythiaInput, noonOnPythiaAxis}, true);
    pythiaInput->Draw("HIST");
    noonOnPythiaAxis->Draw("HIST SAME");
    TLegend *fluxLegend = MakeLegend(0.46, 0.70, 0.92, 0.87);
    fluxLegend->AddEntry(
        pythiaInput, "PYTHIA generated input", "l");
    fluxLegend->AddEntry(
        noonOnPythiaAxis, "n00n flux projection", "l");
    fluxLegend->Draw();
    AddAnalysisLabel(
        0.18, 0.88,
        Form("Common support: %.1f<E_{#gamma}<%.1f MeV",
             commonMinimum, commonMaximum));

    fluxComparisonCanvas->cd(2);
    SetPadStyle(true, true);
    StyleHistogram(pythiaOverNoon, kMagenta + 2, 20);
    pythiaOverNoon->GetXaxis()->SetRangeUser(
        commonMinimum, commonMaximum);
    pythiaOverNoon->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    pythiaOverNoon->GetYaxis()->SetTitle(
        "PYTHIA input / n00n projection");
    pythiaOverNoon->SetMinimum(std::max(
        1.0e-4, 0.5 * PositiveMinimum({pythiaOverNoon})));
    pythiaOverNoon->SetMaximum(std::max(
        2.0, 2.0 * pythiaOverNoon->GetMaximum()));
    pythiaOverNoon->Draw("HIST");
    TLine *unity = new TLine(
        commonMinimum, 1.0, commonMaximum, 1.0);
    unity->SetLineColor(kGray + 2);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->Draw("SAME");
    AddAnalysisLabel(0.18, 0.88, "Shape ratio on common support");
    SaveCanvas(fluxComparisonCanvas, outputPath,
               "13_noon_vs_pythia_egamma");
    extraCanvases.push_back(fluxComparisonCanvas);

    // ------------------------------------------------------------------
    // Persist the additional QA products and canvases.
    // ------------------------------------------------------------------
    const TString outputRootPath =
        outputPath + "/PYTHIA_detector_response_essentials.root";
    TFile outputRoot(outputRootPath, "RECREATE");
    for (TH1 *histogram : outputHistograms) histogram->Write();
    for (TH2 *matrix : outputMatrices) matrix->Write();
    const char *essentialCanvasNames[8] = {
        "c_pythia_weights", "c_pythia_response_matrices",
        "c_pythia_multiplicity_response", "c_pythia_egamma_classes",
        "c_pythia_wgammap", "c_pythia_final_state_qa",
        "c_pythia_correlation_maps", "c_pythia_long_range_dphi"};
    for (const char *name : essentialCanvasNames) {
      TCanvas *canvas = dynamic_cast<TCanvas *>(gROOT->FindObject(name));
      if (canvas) canvas->Write();
    }
    for (TCanvas *canvas : extraCanvases) canvas->Write();
    TNamed(
        "inputFiles",
        Form("PYTHIA=%s; n00n=%s", inputPath.Data(), fluxPath.Data()))
        .Write();
    TNamed(
        "noonProjectionConvention",
        "ProjectionX over all regular impact-parameter bins; interpreted as "
        "stored E_gamma density and normalized over common PYTHIA support")
        .Write();
    outputRoot.Close();

    const TString responsePdf =
        outputPath + "/PYTHIA_detector_response_plots.pdf";
    for (std::size_t index = 0; index < extraCanvases.size(); ++index) {
      TString destination = responsePdf;
      if (index == 0) destination += "(";
      if (index + 1 == extraCanvases.size()) destination += ")";
      extraCanvases[index]->Print(destination);
    }

    std::cout << "\nExtended detector-response plots: "
              << outputPath << std::endl;
    std::cout << "Additional combined PDF: "
              << responsePdf << std::endl;
    std::cout << "Additional ROOT output: "
              << outputRootPath << std::endl;
    std::cout << "n00n/PYTHIA common support: "
              << commonMinimum << " to " << commonMaximum
              << " MeV" << std::endl;

  } catch (const std::exception &error) {
    std::cerr
        << "plot_pythia_post_generation_essentials_macro2 failed: "
        << error.what() << std::endl;
  }

  fluxInput->Close();
  input->Close();
}
