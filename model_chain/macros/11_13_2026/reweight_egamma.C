#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

// Integrate a histogram whose bin contents are densities in E_gamma. Partial
// edge bins are handled using their overlap with [xmin, xmax].
double IntegrateDensity(const TH1 *histogram, double xmin, double xmax)
{
    if (!histogram || !(xmax > xmin)) return 0.0;

    double integral = 0.0;
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
        const double binLow = histogram->GetXaxis()->GetBinLowEdge(bin);
        const double binHigh = histogram->GetXaxis()->GetBinUpEdge(bin);
        const double overlap =
            std::max(0.0, std::min(binHigh, xmax) - std::max(binLow, xmin));
        if (overlap > 0.0) {
            integral += histogram->GetBinContent(bin) * overlap;
        }
    }
    return integral;
}

// Average a density histogram over one energy interval.
void AverageDensityInRange(const TH1 *histogram,
                           double xmin,
                           double xmax,
                           double &average,
                           double &error)
{
    average = 0.0;
    error = 0.0;
    if (!histogram || !(xmax > xmin)) return;

    double integral = 0.0;
    double error2 = 0.0;
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
        const double binLow = histogram->GetXaxis()->GetBinLowEdge(bin);
        const double binHigh = histogram->GetXaxis()->GetBinUpEdge(bin);
        const double overlap =
            std::max(0.0, std::min(binHigh, xmax) - std::max(binLow, xmin));
        if (overlap <= 0.0) continue;

        integral += histogram->GetBinContent(bin) * overlap;
        const double contributionError = histogram->GetBinError(bin) * overlap;
        error2 += contributionError * contributionError;
    }

    const double width = xmax - xmin;
    average = integral / width;
    error = std::sqrt(error2) / width;
}

// Put the target density on the exact PYTHIA energy axis. This avoids
// dividing bins merely because two histograms happen to have the same count.
TH1D *MapDensityToAxis(const TH1 *targetDensity,
                       const TH1 *axisTemplate,
                       const char *name)
{
    if (!targetDensity || !axisTemplate) return nullptr;

    const TAxis *axis = axisTemplate->GetXaxis();
    TH1D *mapped = nullptr;
    if (axis->GetXbins()->GetSize() > 0) {
        mapped = new TH1D(name, "", axis->GetNbins(), axis->GetXbins()->GetArray());
    } else {
        mapped = new TH1D(name, "", axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
    }
    mapped->SetDirectory(nullptr);
    mapped->Sumw2();

    for (int bin = 1; bin <= mapped->GetNbinsX(); ++bin) {
        const double binLow = mapped->GetXaxis()->GetBinLowEdge(bin);
        const double binHigh = mapped->GetXaxis()->GetBinUpEdge(bin);
        double average = 0.0;
        double error = 0.0;
        AverageDensityInRange(targetDensity, binLow, binHigh, average, error);
        mapped->SetBinContent(bin, average);
        mapped->SetBinError(bin, error);
    }
    return mapped;
}

void ZeroOutsideRange(TH1 *histogram, int firstBin, int lastBin)
{
    if (!histogram) return;
    for (int bin = 0; bin <= histogram->GetNbinsX() + 1; ++bin) {
        if (bin < firstBin || bin > lastBin) {
            histogram->SetBinContent(bin, 0.0);
            histogram->SetBinError(bin, 0.0);
        }
    }
}

double MinimumPositive(const TH1 *first, const TH1 *second, const TH1 *third,
                       int firstBin, int lastBin)
{
    double minimum = std::numeric_limits<double>::max();
    const TH1 *histograms[3] = {first, second, third};
    for (const TH1 *histogram : histograms) {
        if (!histogram) continue;
        for (int bin = firstBin; bin <= lastBin; ++bin) {
            const double value = histogram->GetBinContent(bin);
            if (std::isfinite(value) && value > 0.0) {
                minimum = std::min(minimum, value);
            }
        }
    }
    return minimum == std::numeric_limits<double>::max() ? 1.0e-12 : minimum;
}

double MaximumContent(const TH1 *first, const TH1 *second, const TH1 *third,
                      int firstBin, int lastBin)
{
    double maximum = 0.0;
    const TH1 *histograms[3] = {first, second, third};
    for (const TH1 *histogram : histograms) {
        if (!histogram) continue;
        for (int bin = firstBin; bin <= lastBin; ++bin) {
            maximum = std::max(maximum, histogram->GetBinContent(bin));
        }
    }
    return maximum;
}

} // namespace

void reweight_egamma(
    const char *sourcePath =
        "PYTHIA_root_files/egamma_less1p1gev_ep10x100_test17-18_updated_analysis_macro.root",
    const char *targetPath = "n00n_flux_from_Ashik_da/noon_flux.root",
    const char *weightPath =
        "reweighting_factors_ep10x100_test17-18_updated_analysis_macro.txt",
    double minimumUsableSourceContent = 1.0)
{
    gStyle->SetOptStat(0);

    TFile *sourceFile = TFile::Open(sourcePath, "READ");
    if (!sourceFile || sourceFile->IsZombie()) {
        std::cerr << "reweight_egamma failed: could not open PYTHIA source file '"
                  << sourcePath << "'.\n";
        return;
    }
    TH1D *sourceStored = dynamic_cast<TH1D *>(sourceFile->Get("egamma"));
    if (!sourceStored) {
        std::cerr << "reweight_egamma failed: histogram 'egamma' is missing from '"
                  << sourcePath << "'.\n";
        sourceFile->Close();
        return;
    }
    TH1D *sourceCounts = dynamic_cast<TH1D *>(sourceStored->Clone("pythia_source_counts"));
    sourceCounts->SetDirectory(nullptr);
    sourceFile->Close();

    TFile *targetFile = TFile::Open(targetPath, "READ");
    if (!targetFile || targetFile->IsZombie()) {
        std::cerr << "reweight_egamma failed: could not open n00n target file '"
                  << targetPath << "'.\n";
        delete sourceCounts;
        return;
    }
    TH2D *target2D = dynamic_cast<TH2D *>(targetFile->Get("h2fluxNoon_wPnohad"));
    if (!target2D) {
        std::cerr << "reweight_egamma failed: histogram 'h2fluxNoon_wPnohad' "
                  << "is missing from '" << targetPath << "'.\n";
        targetFile->Close();
        delete sourceCounts;
        return;
    }

    // The stored n00n quantity is treated as a density in E_gamma. ProjectionX
    // integrates over the stored regular impact-parameter bins, as in the
    // original macro. Its physical b-integration convention must be validated
    // separately against the producer of noon_flux.root.
    TH1D *targetProjection = target2D->ProjectionX(
        "noon_target_projection", 1, target2D->GetNbinsY(), "e");
    targetProjection->SetDirectory(nullptr);
    targetFile->Close();

    const double targetAxisMin = targetProjection->GetXaxis()->GetXmin();
    const double targetAxisMax = targetProjection->GetXaxis()->GetXmax();

    // Find first and last PYTHIA bins with usable statistics that overlap the
    // n00n energy axis. Their edges define the common support.
    int firstCommonBin = -1;
    int lastCommonBin = -1;
    for (int bin = 1; bin <= sourceCounts->GetNbinsX(); ++bin) {
        const double binLow = sourceCounts->GetXaxis()->GetBinLowEdge(bin);
        const double binHigh = sourceCounts->GetXaxis()->GetBinUpEdge(bin);
        const double content = sourceCounts->GetBinContent(bin);
        const bool overlapsTarget = binHigh > targetAxisMin && binLow < targetAxisMax;
        const bool usable = std::isfinite(content) &&
                            content >= minimumUsableSourceContent;
        if (overlapsTarget && usable) {
            if (firstCommonBin < 0) firstCommonBin = bin;
            lastCommonBin = bin;
        }
    }

    if (firstCommonBin < 0 || lastCommonBin < firstCommonBin) {
        std::cerr << "reweight_egamma failed: no PYTHIA bins with at least "
                  << minimumUsableSourceContent
                  << " entries overlap the n00n target axis.\n";
        delete sourceCounts;
        delete targetProjection;
        return;
    }

    const double commonMinimum = std::max(
        sourceCounts->GetXaxis()->GetBinLowEdge(firstCommonBin), targetAxisMin);
    const double commonMaximum = std::min(
        sourceCounts->GetXaxis()->GetBinUpEdge(lastCommonBin), targetAxisMax);

    TH1D *targetDensity = MapDensityToAxis(
        targetProjection, sourceCounts, "noon_target_normalized_density");
    TH1D *sourceDensity = dynamic_cast<TH1D *>(
        sourceCounts->Clone("pythia_source_normalized_density"));
    sourceDensity->SetDirectory(nullptr);

    ZeroOutsideRange(sourceDensity, firstCommonBin, lastCommonBin);
    ZeroOutsideRange(targetDensity, firstCommonBin, lastCommonBin);

    // PYTHIA stores counts per bin. Divide by bin width after unit-area
    // normalization to obtain a probability density per MeV.
    const double sourceCommonCounts =
        sourceDensity->Integral(firstCommonBin, lastCommonBin);
    // targetDensity already stores density per MeV, so integrate with bin width.
    const double targetCommonIntegral =
        targetDensity->Integral(firstCommonBin, lastCommonBin, "width");

    if (!(sourceCommonCounts > 0.0) || !(targetCommonIntegral > 0.0)) {
        std::cerr << "reweight_egamma failed: common-support integral is non-positive.\n";
        delete sourceCounts;
        delete targetProjection;
        delete sourceDensity;
        delete targetDensity;
        return;
    }

    sourceDensity->Scale(1.0 / sourceCommonCounts, "width");
    targetDensity->Scale(1.0 / targetCommonIntegral);

    // Report the part of the complete n00n spectrum that PYTHIA cannot sample
    // because it lies below the first usable source bin.
    const double targetFullIntegral =
        IntegrateDensity(targetProjection, targetAxisMin, targetAxisMax);
    const double targetBelowThreshold =
        IntegrateDensity(targetProjection, targetAxisMin, commonMinimum);
    const double fractionBelowThreshold =
        targetFullIntegral > 0.0 ? targetBelowThreshold / targetFullIntegral : 0.0;

    std::cout << std::setprecision(8)
              << "\n--- E_gamma common-support report ---\n"
              << "Minimum usable PYTHIA bin content: "
              << minimumUsableSourceContent << "\n"
              << "First usable PYTHIA bin: " << firstCommonBin << "\n"
              << "Last usable PYTHIA bin:  " << lastCommonBin << "\n"
              << "Common support: " << commonMinimum << " < E_gamma < "
              << commonMaximum << " MeV\n"
              << "Fraction of the full n00n projection below the PYTHIA "
              << "threshold: " << fractionBelowThreshold
              << " (" << 100.0 * fractionBelowThreshold << "%)\n"
              << "---------------------------------------\n";

    TH1D *weights = dynamic_cast<TH1D *>(
        targetDensity->Clone("egamma_bin_by_bin_weights"));
    weights->Reset("ICES");
    weights->SetDirectory(nullptr);
    weights->SetTitle("n00n/PYTHIA bin-by-bin weights;E_{#gamma} (MeV);Weight");

    int zeroSourceBinsInsideSupport = 0;
    double uncoveredTargetProbability = 0.0;
    for (int bin = firstCommonBin; bin <= lastCommonBin; ++bin) {
        const double source = sourceDensity->GetBinContent(bin);
        const double target = targetDensity->GetBinContent(bin);
        if (source > 0.0 && std::isfinite(source) && std::isfinite(target)) {
            weights->SetBinContent(bin, target / source);
        } else {
            weights->SetBinContent(bin, 0.0);
            if (target > 0.0) {
                ++zeroSourceBinsInsideSupport;
                uncoveredTargetProbability +=
                    target * targetDensity->GetXaxis()->GetBinWidth(bin);
            }
        }
        weights->SetBinError(bin, 0.0);
    }

    std::cout << "Zero-source bins inside the common range: "
              << zeroSourceBinsInsideSupport << "\n"
              << "Normalized n00n probability in those uncovered bins: "
              << uncoveredTargetProbability << "\n";

    std::ofstream weightFile(weightPath);
    if (!weightFile.is_open()) {
        std::cerr << "reweight_egamma failed: could not write '"
                  << weightPath << "'.\n";
        delete sourceCounts;
        delete targetProjection;
        delete sourceDensity;
        delete targetDensity;
        delete weights;
        return;
    }
    weightFile << "Bin\tEnergy_MeV\tWeight\n";
    weightFile << std::setprecision(12);
    for (int bin = 1; bin <= weights->GetNbinsX(); ++bin) {
        weightFile << bin << '\t'
                   << weights->GetXaxis()->GetBinCenter(bin) << '\t'
                   << weights->GetBinContent(bin) << '\n';
    }
    weightFile.close();
    std::cout << "Weights written to " << weightPath << "\n";

    // Closure: apply each weight to the original PYTHIA count in that bin.
    TH1D *reweightedCounts = dynamic_cast<TH1D *>(
        sourceCounts->Clone("pythia_reweighted_counts"));
    reweightedCounts->SetDirectory(nullptr);
    reweightedCounts->Reset("ICES");
    for (int bin = firstCommonBin; bin <= lastCommonBin; ++bin) {
        const double weight = weights->GetBinContent(bin);
        reweightedCounts->SetBinContent(
            bin, sourceCounts->GetBinContent(bin) * weight);
        reweightedCounts->SetBinError(
            bin, sourceCounts->GetBinError(bin) * std::abs(weight));
    }

    TH1D *reweightedDensity = dynamic_cast<TH1D *>(
        reweightedCounts->Clone("pythia_reweighted_normalized_density"));
    reweightedDensity->SetDirectory(nullptr);
    const double reweightedIntegral =
        reweightedDensity->Integral(firstCommonBin, lastCommonBin);
    if (!(reweightedIntegral > 0.0)) {
        std::cerr << "reweight_egamma failed: reweighted closure integral is non-positive.\n";
        return;
    }
    reweightedDensity->Scale(1.0 / reweightedIntegral, "width");

    double maximumRelativeClosureDifference = 0.0;
    for (int bin = firstCommonBin; bin <= lastCommonBin; ++bin) {
        const double target = targetDensity->GetBinContent(bin);
        if (target <= 0.0 || sourceDensity->GetBinContent(bin) <= 0.0) continue;
        const double difference = std::abs(
            reweightedDensity->GetBinContent(bin) - target) / target;
        maximumRelativeClosureDifference =
            std::max(maximumRelativeClosureDifference, difference);
    }
    std::cout << "Maximum relative closure difference in covered bins: "
              << maximumRelativeClosureDifference << "\n";

    sourceDensity->SetTitle(
        "Before and after n00n reweighting;E_{#gamma} (MeV);"
        "Unit-normalized density (MeV^{-1})");
    sourceDensity->SetLineColor(kBlue + 1);
    sourceDensity->SetLineWidth(2);
    sourceDensity->SetLineStyle(2);
    targetDensity->SetLineColor(kBlack);
    targetDensity->SetLineWidth(3);
    reweightedDensity->SetLineColor(kRed + 1);
    reweightedDensity->SetLineWidth(2);

    const double minimum = MinimumPositive(
        sourceDensity, reweightedDensity, targetDensity,
        firstCommonBin, lastCommonBin);
    const double maximum = MaximumContent(
        sourceDensity, reweightedDensity, targetDensity,
        firstCommonBin, lastCommonBin);
    sourceDensity->GetXaxis()->SetRangeUser(commonMinimum, commonMaximum);
    sourceDensity->SetMinimum(std::max(1.0e-14, 0.3 * minimum));
    sourceDensity->SetMaximum(5.0 * maximum);

    TCanvas *canvas = new TCanvas(
        "c_reweight_egamma_closure", "E_gamma reweighting closure", 1350, 600);
    canvas->Divide(2, 1);

    canvas->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.13);
    sourceDensity->Draw("HIST");
    targetDensity->Draw("HIST SAME");
    reweightedDensity->Draw("HIST SAME");
    TLegend *spectrumLegend = new TLegend(0.46, 0.69, 0.89, 0.88);
    spectrumLegend->SetBorderSize(0);
    spectrumLegend->SetFillStyle(0);
    spectrumLegend->AddEntry(sourceDensity, "PYTHIA before weighting", "l");
    spectrumLegend->AddEntry(reweightedDensity, "PYTHIA after weighting", "l");
    spectrumLegend->AddEntry(targetDensity, "n00n target projection", "l");
    spectrumLegend->Draw();

    TH1D *beforeRatio = dynamic_cast<TH1D *>(
        sourceDensity->Clone("pythia_before_over_noon"));
    TH1D *afterRatio = dynamic_cast<TH1D *>(
        reweightedDensity->Clone("pythia_after_over_noon"));
    beforeRatio->SetDirectory(nullptr);
    afterRatio->SetDirectory(nullptr);
    beforeRatio->Reset("ICES");
    afterRatio->Reset("ICES");

    double ratioMinimum = std::numeric_limits<double>::max();
    double ratioMaximum = 0.0;
    for (int bin = firstCommonBin; bin <= lastCommonBin; ++bin) {
        const double target = targetDensity->GetBinContent(bin);
        if (target <= 0.0) continue;
        const double before = sourceDensity->GetBinContent(bin) / target;
        const double after = reweightedDensity->GetBinContent(bin) / target;
        beforeRatio->SetBinContent(bin, before);
        afterRatio->SetBinContent(bin, after);
        if (before > 0.0) {
            ratioMinimum = std::min(ratioMinimum, before);
            ratioMaximum = std::max(ratioMaximum, before);
        }
        if (after > 0.0) {
            ratioMinimum = std::min(ratioMinimum, after);
            ratioMaximum = std::max(ratioMaximum, after);
        }
    }
    if (ratioMinimum == std::numeric_limits<double>::max()) ratioMinimum = 0.1;
    if (!(ratioMaximum > 0.0)) ratioMaximum = 10.0;

    beforeRatio->SetTitle(
        "Closure ratios;E_{#gamma} (MeV);Spectrum / n00n target");
    beforeRatio->SetLineColor(kBlue + 1);
    beforeRatio->SetLineWidth(2);
    beforeRatio->SetLineStyle(2);
    afterRatio->SetLineColor(kRed + 1);
    afterRatio->SetLineWidth(2);
    beforeRatio->GetXaxis()->SetRangeUser(commonMinimum, commonMaximum);
    beforeRatio->SetMinimum(std::max(1.0e-4, 0.5 * ratioMinimum));
    beforeRatio->SetMaximum(std::max(2.0, 2.0 * ratioMaximum));

    canvas->cd(2);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.13);
    beforeRatio->Draw("HIST");
    afterRatio->Draw("HIST SAME");
    TLine *unity = new TLine(commonMinimum, 1.0, commonMaximum, 1.0);
    unity->SetLineColor(kGray + 2);
    unity->SetLineStyle(3);
    unity->SetLineWidth(2);
    unity->Draw("SAME");
    TLegend *ratioLegend = new TLegend(0.49, 0.73, 0.89, 0.88);
    ratioLegend->SetBorderSize(0);
    ratioLegend->SetFillStyle(0);
    ratioLegend->AddEntry(beforeRatio, "Before / n00n", "l");
    ratioLegend->AddEntry(afterRatio, "After / n00n", "l");
    ratioLegend->AddEntry(unity, "Perfect closure", "l");
    ratioLegend->Draw();

    canvas->SaveAs("reweight_egamma_ep10x100_test17-18_updated_analysis_macro.pdf");
    canvas->SaveAs("reweight_egamma_ep10x100_test17-18_updated_analysis_macro.png");

    TFile outputFile(
        "reweight_egamma_ep10x100_test17-18_updated_analysis_macro.root",
        "RECREATE");
    sourceCounts->Write();
    targetProjection->Write();
    sourceDensity->Write();
    targetDensity->Write();
    weights->Write();
    reweightedCounts->Write();
    reweightedDensity->Write();
    beforeRatio->Write();
    afterRatio->Write();
    canvas->Write();
    outputFile.Close();

    std::cout
        << "Closure plots written to "
        << "reweight_egamma_ep10x100_test17-18_updated_analysis_macro.pdf/.png\n"
        << "Diagnostic histograms written to "
        << "reweight_egamma_ep10x100_test17-18_updated_analysis_macro.root\n";
}
