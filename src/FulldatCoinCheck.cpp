// Function to treat coincidences between D2B and D1B detectors
// Author N. Franel

// ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLine.h>
#include <Rtypes.h>
#include <TApplication.h>
#include "CorrectionTimes.h"
#include <TChain.h>

// C++ headers
#include <iostream>
#include <vector>
#include <stdexcept>
//#include <ROOT/RDataFrame.hxx>
#include <string>


//void FulldatCoinCheck(const std::string& fname_maud = "563_20240623T010937Z_D2B_20240623_195132_0002.sum.root",
//                      const std::string& fname_dssd = "563_20240623T010937Z_D1B_20240623_195133_0002.sum.root")
int main()
{
//    const std::string &fname_maud = "Kiruna_data/563_20240623T010937Z_D2B_20240623_195132_0002.sum.root";
//    const std::string &fname_dssd = "Kiruna_data/563_20240623T010937Z_D1B_20240623_195133_0002.sum.root";
    const std::string &fname_maud = "Kiruna_data/maud_data/*.root";
    const std::string &fname_dssd = "Kiruna_data/dssd_data/*.root";

    std::string fname_maud_corr = "corr_maud_file.root";
    std::string fname_dssd_corr = "corr_dssd_file.root";

    std::string fname_maud_corr_abs = "corr_abs_maud_file.root";
    std::string fname_dssd_corr_abs = "corr_abs_dssd_file.root";

// ====================================================================================================
// Correcting Maud times
// ====================================================================================================
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Opening and correcting Maud times "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
    uint32_t glitch_corr_count_maud = 0;
    uint32_t minor_corr_count_maud = 0;
    CorrectTimes(fname_maud, fname_maud_corr, "maud", glitch_corr_count_maud, minor_corr_count_maud);
    std::cout << " Correction done ! "<< std::endl;

// ====================================================================================================
// Correcting DSSD times
// ====================================================================================================
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Opening and correcting DSSD times "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
    uint32_t glitch_corr_count_dssd = 0;
    uint32_t minor_corr_count_dssd = 0;
    CorrectTimes(fname_dssd, fname_dssd_corr, "dssd", glitch_corr_count_dssd, minor_corr_count_dssd);
    std::cout << " Correction done ! "<< std::endl;

// ====================================================================================================
// creating new timestamps with absolute time reference being pps_info = 1719079200 (22nd june 20:00:00)
// ====================================================================================================
    uint32_t abs_gps_time_ref = 1719079200;
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Absolute time origin is set the 23rd of june 20:00:00  :   " << abs_gps_time_ref << std::endl;
    std::cout << " Time correction to the GPS time is :   " << abs_gps_time_ref << std::endl;
    std::cout << " This value was chosen so that no event has a GPS or PPS time < 0" << std::endl;
    std::cout << "=======================================================================================" << std::endl;
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Assigning the same time origin to PPS is done by setting PPS = GPS at the first entry where GPS time is 10852" << std::endl;
    std::cout << " At this time all instruments are supposed to be working, PPS correction is given and should be the same for all instruments" << std::endl;

    std::cout << "=======================================================================================" << std::endl;
    std::cout << "       Aligning Maud with time origin "<< std::endl;
    ChangeTimeOrigin(fname_maud_corr, fname_maud_corr_abs, "maud", abs_gps_time_ref);
    std::cout << " Time origin set ! "<< std::endl;

    std::cout << "=======================================================================================" << std::endl;
    std::cout << "       Aligning DSSD with time origin "<< std::endl;
    ChangeTimeOrigin(fname_dssd_corr, fname_dssd_corr_abs, "dssd", abs_gps_time_ref);
    std::cout << " Time origin set ! "<< std::endl;

// ====================================================================================================
// Extracting the final root trees
// ====================================================================================================
    auto final_file_maud = new TFile(fname_maud_corr_abs.c_str());
    if (!final_file_maud || final_file_maud->IsZombie()) {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
    auto final_tree_maud = (TTree*) final_file_maud->Get("Events");
//    final_tree_maud->Show(0);

    auto final_file_dssd = new TFile(fname_dssd_corr_abs.c_str());
    if (!final_file_dssd || final_file_dssd->IsZombie()) {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
    auto final_tree_dssd = (TTree*) final_file_dssd->Get("Events");
//    final_tree_dssd->Show(0);

// ====================================================================================================
// Creating tree variables for maud
// ====================================================================================================
    uint32_t ts_init_maud;
    uint32_t pps_cpt_init_maud;
    uint32_t pps_cpt_corr_maud;
    uint32_t gps_corr_maud;
    uint32_t ts_corr_abs_maud;
    uint32_t pps_cpt_corr_abs_maud;
    uint32_t pps_cpt_nocorr_abs_maud;
    uint32_t gps_corr_abs_maud;
    Double_t time_corr_abs_maud;
    final_tree_maud->SetBranchAddress("ts_init",       &ts_init_maud);
    final_tree_maud->SetBranchAddress("pps_cpt_init",  &pps_cpt_init_maud);
    final_tree_maud->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr_maud);
    final_tree_maud->SetBranchAddress("gps_corr",      &gps_corr_maud);
    final_tree_maud->SetBranchAddress("ts_corr_abs",     &ts_corr_abs_maud);
    final_tree_maud->SetBranchAddress("pps_cpt_corr_abs",  &pps_cpt_corr_abs_maud);
    final_tree_maud->SetBranchAddress("pps_cpt_init_abs",  &pps_cpt_nocorr_abs_maud);
    final_tree_maud->SetBranchAddress("gps_corr_abs",      &gps_corr_abs_maud);
    final_tree_maud->SetBranchAddress("time_corr_abs",     &time_corr_abs_maud);

// ====================================================================================================
// Creating tree variables for dssd
// ====================================================================================================
    uint32_t ts_init_dssd;
    uint32_t pps_cpt_init_dssd;
    uint32_t pps_cpt_corr_dssd;
    uint32_t gps_corr_dssd;
    uint32_t ts_corr_abs_dssd;
    uint32_t pps_cpt_corr_abs_dssd;
    uint32_t pps_cpt_nocorr_abs_dssd;
    uint32_t gps_corr_abs_dssd;
    Double_t time_corr_abs_dssd;
    final_tree_dssd->SetBranchAddress("ts_init",       &ts_init_dssd);
    final_tree_dssd->SetBranchAddress("pps_cpt_init",  &pps_cpt_init_dssd);
    final_tree_dssd->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr_dssd);
    final_tree_dssd->SetBranchAddress("gps_corr",      &gps_corr_dssd);
    final_tree_dssd->SetBranchAddress("ts_corr_abs",     &ts_corr_abs_dssd);
    final_tree_dssd->SetBranchAddress("pps_cpt_corr_abs",  &pps_cpt_corr_abs_dssd);
    final_tree_dssd->SetBranchAddress("pps_cpt_init_abs",  &pps_cpt_nocorr_abs_dssd);
    final_tree_dssd->SetBranchAddress("gps_corr_abs",      &gps_corr_abs_dssd);
    final_tree_dssd->SetBranchAddress("time_corr_abs",     &time_corr_abs_dssd);

// ====================================================================================================
// Extracting the number of entries
// ====================================================================================================
    uint64_t nentries_maud = final_tree_maud->GetEntries();
    uint64_t nentries_dssd = final_tree_dssd->GetEntries();
    std::cout << "N events Maud : " << nentries_maud << std::endl;
    std::cout << "N events dssd : " << nentries_dssd << std::endl;

// ====================================================================================================
// Printing the correction statistics
// ====================================================================================================
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Correction applied on Maud timestamps "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
    std::cout << "   " << glitch_corr_count_maud << " over " << nentries_maud << " events : " << static_cast<Double_t>(glitch_corr_count_maud) /static_cast<Double_t>(nentries_maud) * 100. << "% of error" << std::endl;
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Correction applied on DSSD timestamps "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
    std::cout << "   " << glitch_corr_count_dssd << " over " << nentries_dssd << " events : " << static_cast<Double_t>(glitch_corr_count_dssd) /static_cast<Double_t>(nentries_dssd) * 100. << "% of error" << std::endl;


    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Checking the time difference between pps and GPS and comparing Maud and DSSD differences "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
//    uint64_t entry_offset_for_diff = 10; // It's used to account for some misalignment (check is some problems are still there)
    uint32_t old_count;
    uint32_t old_pps_count;
    std::vector<Double_t> vec_gps_maud;
    std::vector<Double_t> diff_gps_pps_maud;
    std::vector<Double_t> diff_gps_pps_maud_corr;
    final_tree_maud->GetEntry(0);
    old_count = gps_corr_abs_maud;
    old_pps_count = pps_cpt_init_maud;
    for (uint64_t i = 0; i<nentries_maud; i++) {
        final_tree_maud->GetEntry(i);
        if (gps_corr_abs_maud != old_count && pps_cpt_init_maud - old_pps_count <= 1){
            vec_gps_maud.push_back(gps_corr_abs_maud);
            diff_gps_pps_maud.push_back(static_cast<Double_t>(pps_cpt_nocorr_abs_maud) - static_cast<Double_t>(gps_corr_abs_maud));
            diff_gps_pps_maud_corr.push_back(static_cast<Double_t>(pps_cpt_corr_abs_maud) - static_cast<Double_t>(gps_corr_abs_maud));
        }
        old_count = gps_corr_abs_maud;
        old_pps_count = pps_cpt_init_maud;
    }
    std::vector<Double_t> vec_gps_dssd;
    std::vector<Double_t> diff_gps_pps_dssd;
    std::vector<Double_t> diff_gps_pps_dssd_corr;
    final_tree_dssd->GetEntry(0);
    old_count = gps_corr_abs_dssd;
    old_pps_count = pps_cpt_init_dssd;
    for (uint64_t i = 0; i<nentries_dssd; i++) {
        final_tree_dssd->GetEntry(i);
        if (gps_corr_abs_dssd != old_count && pps_cpt_init_dssd - old_pps_count <= 1){
            vec_gps_dssd.push_back(gps_corr_abs_dssd);
            diff_gps_pps_dssd.push_back(static_cast<Double_t>(pps_cpt_nocorr_abs_dssd) - static_cast<Double_t>(gps_corr_abs_dssd));
            diff_gps_pps_dssd_corr.push_back(static_cast<Double_t>(pps_cpt_corr_abs_dssd) - static_cast<Double_t>(gps_corr_abs_dssd));
        }
        old_count = gps_corr_abs_dssd;
        old_pps_count = pps_cpt_init_dssd;
    }

// Comparing the difference between the correction of DSSD and Maud
    std::vector<Double_t> maud_dssd_pps_diff_corr_abs;
    std::vector<Double_t> maud_dssd_gps_val_corr_abs;
    for (uint64_t i = 0; i< vec_gps_maud.size(); i++) {
        for (uint64_t j = 0; j< vec_gps_dssd.size(); j++) {
            if (vec_gps_maud[i] == vec_gps_dssd[j]) {
                maud_dssd_pps_diff_corr_abs.push_back(static_cast<Double_t>(diff_gps_pps_maud_corr[i]) - static_cast<Double_t>(diff_gps_pps_dssd_corr[j]));
                maud_dssd_gps_val_corr_abs.push_back(vec_gps_maud[i]);
            }
        }
    }
    std::cout << "     DONE "<< std::endl;

// ====================================================================================================
// Extracting times into a vector (makes things faster)
// ====================================================================================================
    std::vector<Double_t> vec_time_maud;
    for (uint64_t i = 0; i<nentries_maud; i++) {
        final_tree_maud->GetEntry(i);
        vec_time_maud.push_back(time_corr_abs_maud);
    }
    std::vector<Double_t> vec_time_dssd;
    for (uint64_t i = 0; i<nentries_dssd; i++) {
        final_tree_dssd->GetEntry(i);
        vec_time_dssd.push_back(time_corr_abs_dssd);
    }

// ====================================================================================================
// Searching for the coincidences
// ====================================================================================================
    bool display_coinc_vs_window = false;
    if (display_coinc_vs_window) {
        std::cout << "\n=======================================================================================" << std::endl;
        std::cout << " Searching for coincidences between Maud and DSSD "<< std::endl;
        std::cout << "=======================================================================================" << std::endl;

        Double_t window = 1e-8;
        Double_t half_win;
        Double_t step = pow(10, 0.1);
        Double_t win_limit = 1e2;
        uint64_t coinc_counter;
        uint64_t j_restart;
        bool first_ev_in_win;
        std::vector<Double_t> coincidence_num;
        std::vector<Double_t> win_list;
        half_win = window / 2;
        coinc_counter = 0;
        j_restart = 0;
        while (window <= win_limit) {
            std::cout << "Window : " << window << std::endl;
            half_win = window / 2;
            coinc_counter = 0;
            j_restart = 0;
            for (uint64_t i = 0; i<nentries_maud; i++) {
    //            if (i % 1000 == 0) {
    //                std::cout << " i == " << i << " == " << std::endl;
    //            }
                first_ev_in_win = true;
    //            final_tree_maud->GetEntry(i);
                for (uint64_t j = j_restart; j<nentries_dssd; j++) {
    //                final_tree_dssd->GetEntry(j);
                    if (InWindow(vec_time_maud[i], half_win, vec_time_dssd[j])) {
                        coinc_counter += 1;
                        if (first_ev_in_win) {
                            first_ev_in_win = false;
                            j_restart = j;
                        }
    //                    std::cout << i << " " << j << std::endl;
                    }
                    if (vec_time_dssd[j] > vec_time_maud[i] + half_win) {
    //                    std::cout << "dssd : " << vec_time_dssd[j] << " lim maud : " << vec_time_maud[i] + half_win << " bool : " << time_exceed << std::endl;
    //                    std::cout << " == " << std::endl;
                        break;
                    }
                }
            }
            coincidence_num.push_back(coinc_counter);
            win_list.push_back(window);
            std::cout << "  Number of coincidences for a window of " << window << " s is : " << coinc_counter << std::endl;
            window = window * step;
            std::cout << " == " << window << " == " << std::endl;
        }
        auto can4 = new TCanvas("can4", "Coinc vs window");
        can4->SetLogx();
        can4->SetLogy();
        can4->cd(1);
        auto grph10 = new TGraph(win_list.size(), &win_list[0], &coincidence_num[0]);
        grph10->SetTitle("Coinc vs window");
        grph10->SetLineColor(kRed);
        grph10->SetMarkerSize(5.);
        grph10->GetXaxis()->SetTitle("Coincidence window (s)");
        grph10->GetYaxis()->SetTitle("Number of coincidence");
        grph10->Draw("AP");

        can4->Update();
    }

// ====================================================================================================
// printing the time delay between coincident events
// ====================================================================================================
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " printing the time delay between coincident events between Maud and DSSD "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;

    Double_t window_delay = 1e-5;
    Double_t half_win_delay;
    uint64_t coinc_counter_delay;
    uint64_t j_restart_delay;
    bool first_ev_in_win_delay;
    std::vector<Double_t> coincidence_delays;

    half_win_delay = window_delay / 2;
    coinc_counter_delay = 0;
    j_restart_delay = 0;
    for (uint64_t i = 0; i<nentries_maud; i++) {
        first_ev_in_win_delay = true;
//        final_tree_maud->GetEntry(i);
        for (uint64_t j = j_restart_delay; j<nentries_dssd; j++) {
//            final_tree_dssd->GetEntry(j);
            if (InWindow(vec_time_maud[i], half_win_delay, vec_time_dssd[j])) {
                coinc_counter_delay += 1;
                coincidence_delays.push_back(vec_time_maud[i] - vec_time_dssd[j]);
                if (first_ev_in_win_delay) {
                    first_ev_in_win_delay = false;
                    j_restart_delay = j;
                }
//                    std::cout << i << " " << j << std::endl;
            }
            if (vec_time_dssd[j] > vec_time_maud[i] + half_win_delay) {
//                    std::cout << "dssd : " << vec_time_dssd[j] << " lim maud : " << time_corr_abs_maud + half_win_delay << " bool : " << time_exceed << std::endl;
//                    std::cout << " == " << std::endl;
                break;
            }
        }
    }
    std::cout << "  Number of coincidences for a window of " << window_delay << " s is : " << coinc_counter_delay << std::endl;

// ====================================================================================================
// Opening the display application
// ====================================================================================================
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Analysis is done - Opening the display app "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;

    TApplication DisplayApp("App", nullptr, nullptr);
//    auto cant = new TCanvas("cant", "Maud t");
//    final_tree_maud->Draw("pps_cpt_corr_abs:gps_corr_abs");
//    cant->Update();
//    auto can11 = new TCanvas("can11", "Maud timestamp glitches1");
//    can11->Divide(2,2);
//    can11->cd(1);
//    final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 5050 && pps_cpt_init < 5060");
//    can11->cd(2);
//    final_tree_maud->Draw("ts_init/40000000.+pps_cpt_init.:Entry$", "pps_cpt_init > 5050 && pps_cpt_init < 5060");
//    can11->cd(3);
//    final_tree_maud->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_init > 5050 && pps_cpt_init < 5060");
//    can11->cd(4);
//    final_tree_maud->Draw("time_corr_abs.:Entry$", "pps_cpt_init > 5050 && pps_cpt_init < 5060");
//
//    can11->Update();


// ====================================================================================================
// Displaying some values for Maud
// ====================================================================================================
    bool full_trees_display = false;
    if (full_trees_display) {
    // Printing timestamps VS PPS count
        auto can_all1 = new TCanvas("can_all", "All - ts vs pps");
        can_all1->Divide(2,2);
        can_all1->cd(1);
        final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init_abs");
        can_all1->cd(2);
        final_tree_dssd->Draw("ts_init/40000000.:pps_cpt_init_abs");
        can_all1->cd(3);
        final_tree_maud->Draw("ts_corr_abs/40000000.:pps_cpt_corr_abs");
        can_all1->cd(4);
        final_tree_dssd->Draw("ts_corr_abs/40000000.:pps_cpt_corr_abs");

        can_all1->Update();

    // Printing timestamps VS PPS count
        auto can_all2 = new TCanvas("can_all", "All - ts vs pps");
        can_all2->Divide(2,2);
        can_all2->cd(1);
        final_tree_maud->Draw("ts_init/40000000.+pps_cpt_init_abs.:Entry$");
        can_all2->cd(2);
        final_tree_dssd->Draw("ts_init/40000000.+pps_cpt_init_abs.:Entry$");
        can_all2->cd(3);
        final_tree_maud->Draw("ts_corr_abs/40000000.+pps_cpt_corr_abs.:Entry$");
        can_all2->cd(4);
        final_tree_dssd->Draw("ts_corr_abs/40000000.+pps_cpt_corr_abs.:Entry$");

        can_all2->Update();
    }

    bool display_ts_corr_maud = false;
    if (display_ts_corr_maud) {
        // Display of the first glitch for 77560 < pps_cpt < 77580
        auto can11 = new TCanvas("can11", "Maud timestamp glitches1");
        can11->Divide(2,2);
        can11->cd(1);
        final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 77560 && pps_cpt_init < 77580");
        can11->cd(2);
        final_tree_maud->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 77560 && pps_cpt_init < 77580");
        can11->cd(3);
        final_tree_maud->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_corr > 77560 && pps_cpt_corr < 77580");
        can11->cd(4);
        final_tree_maud->Draw("time_corr_abs-76780.:Entry$", "pps_cpt_corr > 77560 && pps_cpt_corr < 77580");

        can11->Update();

        // Display of the second glitch for 79770 < pps_cpt < 79850
        auto can12 = new TCanvas("can12", "Maud timestamp glitches2");
        can12->Divide(1,2);
        can12->cd(1);
        final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 79770 && pps_cpt_init < 79850");
        can12->cd(2);
        final_tree_maud->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 79770 && pps_cpt_init < 79850");

        can12->Update();

        // Display of the 3rd glitch for 83490 < pps_cpt < 83510
        auto can13 = new TCanvas("can13", "Maud timestamp glitches3");
        can13->Divide(2,2);
        can13->cd(1);
        final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 83490 && pps_cpt_init < 83510");
        can13->cd(2);
        final_tree_maud->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 83490 && pps_cpt_init < 83510");
        can13->cd(3);
        final_tree_maud->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_corr > 83490 && pps_cpt_corr < 83510");
        can13->cd(4);
        final_tree_maud->Draw("time_corr_abs-76780.:Entry$", "pps_cpt_corr > 83490 && pps_cpt_corr < 83510");

        can13->Update();

        // Display of the 4th glitch for 84170 < pps_cpt < 84180
        auto can14 = new TCanvas("can14", "Maud timestamp glitches4");
        can14->Divide(2,2);
        can14->cd(1);
        final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 84170 && pps_cpt_init < 84180");
        can14->cd(2);
        final_tree_maud->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 84170 && pps_cpt_init < 84180");
        can14->cd(3);
        final_tree_maud->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_corr > 84170 && pps_cpt_corr < 84180");
        can14->cd(4);
        final_tree_maud->Draw("time_corr_abs-76780.:Entry$", "pps_cpt_corr > 84170 && pps_cpt_corr < 84180");

        can14->Update();

        // Display of the whole file
        auto can15 = new TCanvas("can15", "Maud timestamps");
        can15->Divide(2,1);
        can15->cd(1);
        final_tree_maud->Draw("ts_init/40000000.:pps_cpt_init");
        can15->cd(2);
        final_tree_maud->Draw("ts_corr_abs/40000000.:pps_cpt_corr");

        can15->Update();
    }
    bool display_comp_clocks_maud = false;
    if (display_comp_clocks_maud) {
        auto can16 = new TCanvas("can16", "Comp clocks Maud");
        can16->Divide(2,2);

        can16->cd(1);
        final_tree_maud->Draw("pps_cpt_init:gps_corr");

        can16->cd(2);
        final_tree_maud->Draw("pps_cpt_corr:gps_corr");

        can16->cd(3);
        auto grph3 = new TGraph(vec_gps_maud.size(), &vec_gps_maud[0], &diff_gps_pps_maud[0]);
        grph3->SetTitle("diff GPS PPS");
        grph3->GetXaxis()->SetTitle("GPS time (s)");
        grph3->GetYaxis()->SetTitle("time diff (s)");
        grph3->Draw("AP");

        can16->cd(4);
        auto grph4 = new TGraph(vec_gps_maud.size(), &vec_gps_maud[0], &diff_gps_pps_maud_corr[0]);
        grph4->SetTitle("diff GPS pps corr");
        grph4->GetXaxis()->SetTitle("GPS time (s)");
        grph4->GetYaxis()->SetTitle("time diff corrected (s)");
        grph4->Draw("AP");

        can16->Update();

//        auto grph1 = new TGraph(abs_gps_maud.size(), &abs_gps_maud[0], &abs_pps_maud[0]);
//        grph1->SetTitle("GPS vs PPS");
//        grph1->SetLineColor(kRed);
//        grph1->GetXaxis()->SetTitle("GPS time (s)");
//        grph1->GetYaxis()->SetTitle("PPS time (s)");
//        grph1->Draw("AP");
//        auto grph2 = new TGraph(abs_gps_maud.size(), &abs_gps_maud[0], &abs_pps_maud_corrected[0]);
//        grph2->SetTitle("GPS vs PPS corr");
//        grph2->SetLineColor(kBlue);
//        grph2->GetXaxis()->SetTitle("GPS time (s)");
//        grph2->GetYaxis()->SetTitle("PPS corrected time (s)");
//        grph2->Draw("AP");
    }
// ====================================================================================================
// Displaying some values for DSSD
// ====================================================================================================
    bool display_ts_corr_dssd = false;
    if (display_ts_corr_dssd) {
        // Display of the first glitch for 77560 < pps_cpt < 77580
        auto can21 = new TCanvas("can21", "DSSD timestamp glitches1");
        can21->Divide(2,2);
        can21->cd(1);
        final_tree_dssd->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 77560 && pps_cpt_init < 77580");
        can21->cd(2);
        final_tree_dssd->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 77560 && pps_cpt_init < 77580");
        can21->cd(3);
        final_tree_dssd->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_corr > 77560 && pps_cpt_corr < 77580");
        can21->cd(4);
        final_tree_dssd->Draw("time_corr_abs-76780.:Entry$", "pps_cpt_corr > 77560 && pps_cpt_corr < 77580");

        can21->Update();

        // Display of the second glitch for 79770 < pps_cpt < 79850
        auto can22 = new TCanvas("can22", "DSSD timestamp glitches2");
        can22->Divide(1,2);
        can22->cd(1);
        final_tree_dssd->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 79770 && pps_cpt_init < 79850");
        can22->cd(2);
        final_tree_dssd->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 79770 && pps_cpt_init < 79850");

        can22->Update();

        // Display of the 3rd glitch for 83490 < pps_cpt < 83510
        auto can23 = new TCanvas("can23", "DSSD timestamp glitches3");
        can23->Divide(2,2);
        can23->cd(1);
        final_tree_dssd->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 83490 && pps_cpt_init < 83510");
        can23->cd(2);
        final_tree_dssd->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 83490 && pps_cpt_init < 83510");
        can23->cd(3);
        final_tree_dssd->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_corr > 83490 && pps_cpt_corr < 83510");
        can23->cd(4);
        final_tree_dssd->Draw("time_corr_abs-76780.:Entry$", "pps_cpt_corr > 83490 && pps_cpt_corr < 83510");

        can23->Update();

        // Display of the 4th glitch for 84170 < pps_cpt < 84180
        auto can24 = new TCanvas("can24", "DSSD timestamp glitches4");
        can24->Divide(2,2);
        can24->cd(1);
        final_tree_dssd->Draw("ts_init/40000000.:pps_cpt_init", "pps_cpt_init > 84170 && pps_cpt_init < 84180");
        can24->cd(2);
        final_tree_dssd->Draw("ts_init/40000000.+pps_cpt_init-76780.:Entry$", "pps_cpt_init > 84170 && pps_cpt_init < 84180");
        can24->cd(3);
        final_tree_dssd->Draw("ts_corr_abs/40000000.:pps_cpt_corr", "pps_cpt_corr > 84170 && pps_cpt_corr < 84180");
        can24->cd(4);
        final_tree_dssd->Draw("time_corr_abs-76780.:Entry$", "pps_cpt_corr > 84170 && pps_cpt_corr < 84180");

        can24->Update();

        // Display of the whole file
        auto can25 = new TCanvas("can25", "DSSD timestamps");
        can25->Divide(2,1);
        can25->cd(1);
        final_tree_dssd->Draw("ts_init/40000000.:pps_cpt_init");
        can25->cd(2);
        final_tree_dssd->Draw("ts_corr_abs/40000000.:pps_cpt_corr");

        can25->Update();
    }
    bool display_comp_clocks_dssd = false;
    if (display_comp_clocks_dssd) {
        auto can26 = new TCanvas("can26", "Comp clocks DSSD");
        can26->Divide(2,2);

        can26->cd(1);
        final_tree_dssd->Draw("pps_cpt_init:gps_corr");

        can26->cd(2);
        final_tree_dssd->Draw("pps_cpt_corr:gps_corr");

        can26->cd(3);
        auto grph7 = new TGraph(vec_gps_dssd.size(), &vec_gps_dssd[0], &diff_gps_pps_dssd[0]);
        grph7->SetTitle("diff GPS PPS");
        grph7->GetXaxis()->SetTitle("GPS time (s)");
        grph7->GetYaxis()->SetTitle("time diff (s)");
        grph7->Draw("AP");

        can26->cd(4);
        auto grph8 = new TGraph(vec_gps_dssd.size(), &vec_gps_dssd[0], &diff_gps_pps_dssd_corr[0]);
        grph8->SetTitle("diff GPS pps corr");
        grph8->GetXaxis()->SetTitle("GPS time (s)");
        grph8->GetYaxis()->SetTitle("time diff corrected (s)");
        grph8->Draw("AP");

        can26->Update();

//        auto grph5 = new TGraph(abs_gps_dssd.size(), &abs_gps_dssd[0], &abs_pps_dssd[0]);
//        grph5->SetTitle("GPS vs PPS");
//        grph5->SetLineColor(kRed);
//        grph5->GetXaxis()->SetTitle("GPS time (s)");
//        grph5->GetYaxis()->SetTitle("PPS time (s)");
//        grph5->Draw("AP");
//        auto grph6 = new TGraph(abs_gps_dssd.size(), &abs_gps_dssd[0], &abs_pps_dssd_corrected[0]);
//        grph6->SetTitle("GPS vs PPS corr");
//        grph6->SetLineColor(kBlue);
//        grph6->GetXaxis()->SetTitle("GPS time (s)");
//        grph6->GetYaxis()->SetTitle("PPS corrected time (s)");
//        grph6->Draw("AP");
    }
    bool display_comp_clocks_diff = false;
    if (display_comp_clocks_diff) {
        auto can3 = new TCanvas("can3", "Comp of clock diffs DSSD & Maud");
        can3->cd(1);
        auto grph9 = new TGraph(maud_dssd_gps_val_corr_abs.size(), &maud_dssd_gps_val_corr_abs[0], &maud_dssd_pps_diff_corr_abs[0]);
        grph9->SetTitle("Diff of clock diff");
        grph9->SetLineColor(kRed);
        grph9->GetXaxis()->SetTitle("GPS time (s)");
        grph9->GetYaxis()->SetTitle("Time diff (s)");
        grph9->Draw("AP");

        can3->Update();
    }
    bool display_coinc_delay = false;
    if (display_coinc_delay) {
        auto can5 = new TCanvas("can5", "Coinc delay histogram");
//        can4->SetLogx();
        can5->SetLogy();
        int bin_num = half_win_delay / 2e-8;
        std::vector<double> weights(coincidence_delays.size(),1);
        can5->cd(1);
        std::cout << "n bins : " << bin_num << std::endl;
        std::cout << "size coinc list : " << coincidence_delays.size() << std::endl;
        auto h11 = new TH1D("h11", "Coinc delay histogram", bin_num, -half_win_delay,  half_win_delay);
        h11->FillN(coincidence_delays.size(), coincidence_delays.data(), weights.data());
//        h11->SetTitle("Coinc vs window");
//        h11->SetLineColor(kRed);
//        h11->SetMarkerSize(5.);
//        h11->GetXaxis()->SetTitle("Coincidence window (s)");
//        h11->GetYaxis()->SetTitle("Number of coincidence");
        h11->Draw();

        can5->Update();
    }
    DisplayApp.Run();
    return 0;
}
