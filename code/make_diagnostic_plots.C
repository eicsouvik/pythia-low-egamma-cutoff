#include <algorithm>
#include <stdexcept>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

void make_diagnostic_plots(
    const char *inputPath =
        "inputs/egamma_less1p1gev_ep10x100_test17-18_eta2_01.root",
    const char *outputStem = "outputs/low_energy_cutoff_QA") {
  gStyle->SetOptStat(0);

  TFile input(inputPath, "READ");
  if (input.IsZombie()) {
    throw std::runtime_error(Form("Could not open %s", inputPath));
  }

  TH1 *egamma = dynamic_cast<TH1 *>(input.Get("egamma"));
  TH1 *npartRaw = dynamic_cast<TH1 *>(input.Get("h_npart_raw"));
  if (!egamma || !npartRaw) {
    throw std::runtime_error(
        "Input must contain 'egamma' and 'h_npart_raw' histograms");
  }

  TH1 *egammaPlot = dynamic_cast<TH1 *>(egamma->Clone("egamma_diagnostic"));
  TH1 *npartPlot = dynamic_cast<TH1 *>(npartRaw->Clone("npart_raw_diagnostic"));
  egammaPlot->SetDirectory(nullptr);
  npartPlot->SetDirectory(nullptr);

  egammaPlot->SetTitle("Generated photon spectrum;E_{#gamma} (MeV);Events");
  npartPlot->SetTitle("Raw final-state multiplicity;N_{part}^{raw};Events");
  egammaPlot->SetLineColor(kBlue + 1);
  egammaPlot->SetLineWidth(2);
  npartPlot->SetLineColor(kRed + 1);
  npartPlot->SetLineWidth(2);

  TCanvas canvas("c_low_energy_diagnostics", "Low-energy cutoff diagnostics",
                 1500, 650);
  canvas.Divide(2, 1);

  canvas.cd(1);
  //gPad->SetLogy();
  //gPad->SetLeftMargin(0.13);
  //gPad->SetBottomMargin(0.13);
  egammaPlot->GetXaxis()->SetRangeUser(0., 200.);
  const double egammaMaximum = std::max(1.0, egammaPlot->GetMaximum());
  TH1F *egammaFrame = gPad->DrawFrame(0., 0., 200., 1.10 * egammaMaximum);
  egammaFrame->SetTitle("Generated photon spectrum;E_{#gamma} (MeV);Events");
  egammaPlot->Draw("HIST SAME");
  TLine egammaThreshold(62.5, 0., 62.5, 3e5);
  egammaThreshold.SetLineColor(kBlack);
  egammaThreshold.SetLineStyle(2);
  egammaThreshold.SetLineWidth(2);
  egammaThreshold.Draw();
  TLatex photonLabel;
  photonLabel.SetNDC();
  photonLabel.SetTextSize(0.038);
  photonLabel.DrawLatex(0.48, 0.83,
                       "W_{min}=5 GeV #Rightarrow E_{#gamma}#approx62.5 MeV");

  canvas.cd(2);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  npartPlot->GetXaxis()->SetRangeUser(0., 50.);
  npartPlot->Draw("HIST");
  TLine multiplicityThreshold(4., std::max(1.0, npartPlot->GetMinimum(0.0)),
                              4., npartPlot->GetMaximum());
  multiplicityThreshold.SetLineColor(kBlack);
  multiplicityThreshold.SetLineStyle(2);
  multiplicityThreshold.SetLineWidth(2);
  multiplicityThreshold.Draw();
  TLatex multiplicityLabel;
  multiplicityLabel.SetNDC();
  multiplicityLabel.SetTextSize(0.038);
  multiplicityLabel.DrawLatex(0.48, 0.83,
                             "Observed edge near N_{part}^{raw}=4");

  canvas.SaveAs(Form("%s.pdf", outputStem));
  //canvas.SaveAs(Form("%s.png", outputStem));

  delete egammaPlot;
  delete npartPlot;
}
