#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream> 
#include <fstream>

void reweight_egamma() {
    gStyle->SetOptStat(0);

    // 1. Access the source histogram
    TFile *source_file = TFile::Open("egamma_less1p1gev_ep2x100_test2.root", "READ");
    if (!source_file || source_file->IsZombie()) {
        std::cout << "Error: Could not open source file." << std::endl;
        return;
    }
    TH1D *h_source = (TH1D*)source_file->Get("egamma");
    if (!h_source) {
        std::cout << "Error: Could not find 'egamma' histogram in source file." << std::endl;
        return;
    }
    h_source->SetTitle("Source Distribution (egamma)");
    h_source->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    h_source->GetYaxis()->SetTitle("Counts");

    // 2. Access the target histogram and project it
    TFile *target_file = TFile::Open("noon_flux.root", "READ");
    if (!target_file || target_file->IsZombie()) {
        std::cout << "Error: Could not open target file." << std::endl;
        return;
    }
    TH2D *h2_target_raw = (TH2D*)target_file->Get("h2fluxNoon_wPnohad");
    if (!h2_target_raw) {
        std::cout << "Error: Could not find 'h2fluxNoon_wPnohad' histogram in target file." << std::endl;
        return;
    }
    // Project the 2D histogram to its X-axis to get the target shape
    TH1D *h_target = h2_target_raw->ProjectionX("h_target");
    h_target->SetTitle("Target Distribution (noon_flux projection)");
    h_target->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
    h_target->GetYaxis()->SetTitle("Counts");

    // --- The Reweighting Process ---

    // Ensure histograms are compatible for division
    int nbins_source = h_source->GetNbinsX();
    int nbins_target = h_target->GetNbinsX();

    if (nbins_source != nbins_target) {
        std::cout << "Histograms have different binning (" << nbins_source << " vs " << nbins_target 
                  << "). Attempting to rebin..." << std::endl;
                  
        if (nbins_target > nbins_source && nbins_target % nbins_source == 0) {
            // Target has more bins, so we rebin it to match the source
            int rebin_factor = nbins_target / nbins_source;
            h_target->Rebin(rebin_factor);
            std::cout << "Rebinned target by a factor of " << rebin_factor << std::endl;
        } else if (nbins_source > nbins_target && nbins_source % nbins_target == 0) {
            // Source has more bins, so we rebin it to match the target
            int rebin_factor = nbins_source / nbins_target;
            h_source->Rebin(rebin_factor);
            std::cout << "Rebinned source by a factor of " << rebin_factor << std::endl;
        } else {
            // Bin counts are not integer multiples, so we cannot proceed
            std::cout << "Error: Binning is incompatible and cannot be automatically reconciled." << std::endl;
            return;
        }
    }

    // Normalize both histograms to unit area before creating weights
    h_source->Scale(1.0 / h_source->Integral());
    //h_target->Scale(1.0 / h_target->Integral());

    // 3. Create the weight histogram by dividing Target by Source
    TH1D *h_weights = (TH1D*)h_target->Clone("h_weights");
    h_weights->SetTitle("Calculated Weights");
    h_weights->Divide(h_source);

    std::ofstream weight_file("reweighting_factors_ep2x100_test2.txt");
    if (!weight_file.is_open()) {
        std::cerr << "Error: Could not open reweighting_factors.txt for writing." << std::endl;
    } else {
        weight_file << "Bin\tEnergy_MeV\tWeight" << std::endl; // Write header
        for (int i = 1; i <= h_source->GetNbinsX(); ++i) {
            double bin_center_mev = h_weights->GetXaxis()->GetBinCenter(i);
            double weight = h_weights->GetBinContent(i);

            // Write energy (GeV) and weight, tab-separated
            weight_file << i << "\t" <<bin_center_mev << "\t" << weight << std::endl;
        }
        weight_file.close();
        std::cout << "Reweighting factors saved to reweighting_factors.txt" << std::endl;
    }

    // 4. Create the final, reweighted histogram
    TH1D *h_reweighted = (TH1D*)h_source->Clone("h_reweighted");
    h_reweighted->SetTitle("Reweighted Source Distribution");
    h_reweighted->Reset(); // Clear the cloned contents

    // 5. Loop through the original source histogram bins to apply weights
    
    // --- START MODIFICATION ---
    std::cout << "\n--- Calculated Reweighting Factors ---" << std::endl;
    for (int i = 1; i <= h_source->GetNbinsX(); ++i) {
        double original_content = h_source->GetBinContent(i);
        double weight = h_weights->GetBinContent(i);
        
        // Print the weight for the current bin
        double bin_center_energy = h_weights->GetXaxis()->GetBinCenter(i);
        std::cout << "Bin " << i << " (E_gamma ~ " << bin_center_energy << " MeV): Weight = " << weight << std::endl;

        h_reweighted->SetBinContent(i, original_content * weight);
    }
    std::cout << "-------------------------------------\n" << std::endl;
    // --- END MODIFICATION ---
    
    // Normalize the final reweighted histogram for a clean shape comparison
    h_reweighted->Scale(1.0 / h_reweighted->Integral());

    // --- Visualization ---
    TCanvas *c1 = new TCanvas("c1", "Reweighting EGamma", 1200, 500);
    c1->Divide(2, 1);

    // Original distributions
    c1->cd(1);
    h_source->SetLineColor(kBlue);
    h_source->Draw("hist");
    h_target->SetLineColor(kRed);
    h_target->Draw("hist same");
    TLegend *leg1 = new TLegend(0.5, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_source, "Original Source (egamma)", "l");
    leg1->AddEntry(h_target, "Target Shape (noon_flux)", "l");
    leg1->Draw();

    // Reweighted result
    c1->cd(2);
    h_reweighted->SetLineColor(kBlue);
    h_reweighted->Draw("hist");
    h_target->SetLineColor(kRed);
    h_target->Draw("hist same");
    TLegend *leg2 = new TLegend(0.5, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_reweighted, "Reweighted Source", "l");
    leg2->AddEntry(h_target, "Target Shape", "l");
    leg2->Draw();
}