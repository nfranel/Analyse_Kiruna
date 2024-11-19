#include "CorrectionTimes.h"
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <Rtypes.h>
#include <string>
#include <fstream>


bool InWindow(const Double_t mid_interval, const Double_t half_window, const Double_t tested_value) {
// If       mid_interval - half_window/2 <= tested_value <= mid_interval + half_window/2      returns true, else return false
    if (tested_value >= mid_interval - half_window && tested_value <= mid_interval + half_window) {
//        std::cout << fixed << setprecision(8) << mid_interval - half_window << " <= " << tested_value << "<=" << mid_interval + half_window << std::endl;
        return true;
    } else {
        return false;
    }
}


void file_exists(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.good()) {
        std::cerr << "Error: File " << filename << " does not exist or cannot be accessed." << std::endl;
    }
    return;
}


// Implémentation de la fonction CorrectTimes for Maud
void CorrectTimes(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &glitch_corr_count, uint32_t &minor_corr_count) {
    //
    // Function used to correct the timestamps of the root files.
    //
    // Extraction of the root tree
    TChain *tree_init = new TChain("Events");
    tree_init->Add(init_tree_name.c_str());
    tree_init->Show(100);

    // Variables used to obtain tree values
    uint32_t ts;
    uint32_t pps_cpt;
    uint32_t gps;

    // Declaring file_final
    TFile* file_final = nullptr;
    // Declaring tree_final
    TTree* tree_final = nullptr;

    // Connect to the branches and opening new root tree according to the detector
    if (detector == "maud") {
        // connect Maud branches
        tree_init->SetBranchAddress("timestamp",   &ts);
        tree_init->SetBranchAddress("pps_cpt",     &pps_cpt);
        tree_init->SetBranchAddress("pps_info",    &gps);

        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Corrected root file for Maud");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Corrected Maud tree");
    } else if (detector == "dssd") {
        // connect dssd branches
        tree_init->SetBranchAddress("timestamp",   &ts);
        tree_init->SetBranchAddress("pps_cpt",     &pps_cpt);
        tree_init->SetBranchAddress("pps_info",    &gps);

        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Corrected root file for DSSD");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Corrected DSSD tree");
    } else if (detector == "cea") {
        std::cerr << " Error : cea detector not implemented yet " << std::endl;
    } else if (detector == "ucd") {
        std::cerr << " Error : ucd detector not implemented yet " << std::endl;
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud, dssd, cea or ucd" << std::endl;
    }

    // Variables to create the new tree + linking them to the tree
    uint32_t ts_init;
    uint32_t pps_cpt_init;
    uint32_t ts_corr;
    uint32_t pps_cpt_corr;
    uint32_t gps_corr;
    Double_t time_corr;
    tree_final->Branch("ts_init", &ts_init, "ts_init/i");
    tree_final->Branch("pps_cpt_init", &pps_cpt_init, "pps_cpt_init/i");
    tree_final->Branch("ts_corr", &ts_corr, "ts_corr/i");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr, "pps_cpt_corr/i");
    tree_final->Branch("gps_corr", &gps_corr, "gps_corr/i");
    tree_final->Branch("time_corr", &time_corr, "time_corr/D");
//    tree_final->SetBranchAddress("ts_init",       &ts_init);
//    tree_final->SetBranchAddress("pps_cpt_init",  &pps_cpt_init);
//    tree_final->SetBranchAddress("ts_corr",       &ts_corr);
//    tree_final->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr);
//    tree_final->SetBranchAddress("gps_corr",      &gps_corr);
//    tree_final->SetBranchAddress("time_corr",     &time_corr);

    // Variables not saved, used for correcting the time series and creating the new tree
    uint32_t correction_counter = 0;
    uint32_t old_ts=0;
//    uint32_t read_ts;
    uint32_t next_ts;
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Correction of " << detector << " timestamps "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
    uint64_t nentries = tree_init->GetEntries();
    for (uint64_t i = 0; i<nentries; i++) {
        // Loading the tree entry
        tree_init->GetEntry(i);
        // Getting the initial, uncorrected, tree values
        pps_cpt_init = pps_cpt;
        ts_init = ts;
        // Starting the correction if needed
        if (ts >= 40.e6) {
            ts_corr = ts - 40.e6;
            tree_init->GetEntry(i+1);
            next_ts = ts;
            tree_init->GetEntry(i);
            if (old_ts < 40.e6 && next_ts >= 40.e6) {
                correction_counter += 1;
                pps_cpt_corr = pps_cpt + correction_counter;
                std::cout << detector << " correction on entry " << i << ", pps count " << pps_cpt << " to " << pps_cpt_corr << std::endl;
                glitch_corr_count += 1;
            } else if (old_ts < 40.e6 && next_ts < 40.e6){
                std::cout << "ts threshold exceeded but might be due to delay in saving the value, correction applied for this value only, verification is needed to see if there is no shift compared to the GPS clock" << std::endl;
                pps_cpt_corr = pps_cpt + correction_counter + 1;
                std::cout << " Value for ts and pps before local correction : " << ts << " " << pps_cpt + correction_counter << " and after : " << ts_corr <<  " " << pps_cpt_corr << std::endl;
                tree_init->GetEntry(i+1);
                if (static_cast<Double_t>(ts)/40.e6+static_cast<Double_t>(pps_cpt)+static_cast<Double_t>(correction_counter) > static_cast<Double_t>(ts_corr)/40.e6+static_cast<Double_t>(pps_cpt_corr)) {
                    std::cout << "     Successful correction : the next entry time is greater than the corrected one : " << std::endl;
                    std::cout << "     Next entry has the following ts and pps : " << ts << " " << pps_cpt + correction_counter << std::endl;
                } else {
                    std::cerr << "     Correction failed : the next entry has a time smaller than the corrected time - Another correction should be used or suppressing the value. " << std::endl;
                    std::cerr << "     Value of next entry (not corrected) : " << ts << " " << pps_cpt + correction_counter << std::endl;
                }
                tree_init->GetEntry(i);
                minor_corr_count += 1;
            } else {
                pps_cpt_corr = pps_cpt + correction_counter;
            }
        } else if (ts >= 80.e6) {
            throw std::runtime_error("ERROR : 2 pps incrementations not done in a row, code update needed to solve this");
        } else {
            ts_corr = ts;
            pps_cpt_corr = pps_cpt + correction_counter;
        }
        gps_corr = gps;
        time_corr = pps_cpt_corr + ts_corr/40.e6;

        tree_final->Fill();

        old_ts = ts;
    }
//    file_init->Close();
    file_final->Write();
    file_final->Close();
}

void ChangeTimeOrigin(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &abs_gps_time_ref) {
    //
    // Function used to correct the align the timestamp with a common time origin
    //

    // Extraction of the root tree
    file_exists(init_tree_name);
    auto file_init = new TFile(init_tree_name.c_str());
    if (!file_init || file_init->IsZombie()) {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
//    std::cout << "KEYS : " << file_init->ls() << std::endl;
//    file_init->ls();
    auto tree_init = (TTree*) file_init->Get("Events");
//    std::cout << "NOMBRE DENTREES : " << tree_init->GetEntries() << std::endl;

    // Variables to read the initial tree
    uint32_t ts_init;
    uint32_t pps_cpt_init;
    uint32_t ts_corr;
    uint32_t pps_cpt_corr;
    uint32_t gps_corr;
    Double_t time_corr;
    tree_init->SetBranchAddress("ts_init",       &ts_init);
    tree_init->SetBranchAddress("pps_cpt_init",  &pps_cpt_init);
    tree_init->SetBranchAddress("ts_corr",       &ts_corr);
    tree_init->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr);
    tree_init->SetBranchAddress("gps_corr",      &gps_corr);
    tree_init->SetBranchAddress("time_corr",     &time_corr);

    // Declaring file_final
    TFile* file_final = nullptr;
    // Declaring tree_final
    TTree* tree_final = nullptr;

    // Opening new root tree according to the detector
    if (detector == "maud") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for Maud");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time Maud tree");
    } else if (detector == "dssd") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for DSSD");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time DSSD tree");
    } else if (detector == "cea") {
        std::cerr << " Error : cea detector not implemented yet " << std::endl;
    } else if (detector == "ucd") {
        std::cerr << " Error : ucd detector not implemented yet " << std::endl;
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud, dssd, cea or ucd" << std::endl;
    }

    // Variables to create the new tree + linking them to the tree
    uint32_t ts_init_noabs;
    uint32_t pps_cpt_init_noabs;
    uint32_t pps_cpt_corr_noabs;
    uint32_t gps_corr_noabs;
    uint32_t ts_corr_abs;
    uint32_t pps_cpt_corr_abs;
    uint32_t pps_cpt_nocorr_abs;
    uint32_t gps_corr_abs;
    Double_t time_corr_abs;
    tree_final->Branch("ts_init", &ts_init_noabs, "ts_init/i");
    tree_final->Branch("pps_cpt_init", &pps_cpt_init_noabs, "pps_cpt_init/i");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr_noabs, "pps_cpt_corr/i");
    tree_final->Branch("gps_corr", &gps_corr_noabs, "gps_corr/i");
    tree_final->Branch("ts_corr_abs", &ts_corr_abs, "ts_corr_abs/i");
    tree_final->Branch("pps_cpt_corr_abs", &pps_cpt_corr_abs, "pps_cpt_corr_abs/i");
    tree_final->Branch("pps_cpt_init_abs", &pps_cpt_nocorr_abs, "pps_cpt_init_abs/i");
    tree_final->Branch("gps_corr_abs", &gps_corr_abs, "gps_corr_abs/i");
    tree_final->Branch("time_corr_abs", &time_corr_abs, "time_corr_abs/D");

//    tree_final->SetBranchAddress("ts_init",       &ts_init_noabs);
//    tree_final->SetBranchAddress("pps_cpt_init",  &pps_cpt_init_noabs);
//    tree_final->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr_noabs);
//    tree_final->SetBranchAddress("gps_corr",      &gps_corr_noabs);
//    tree_final->SetBranchAddress("ts_corr_abs",     &ts_corr_abs);
//    tree_final->SetBranchAddress("pps_cpt_corr_abs",  &pps_cpt_corr_abs);
//    tree_final->SetBranchAddress("pps_cpt_init_abs",  &pps_cpt_nocorr_abs);
//    tree_final->SetBranchAddress("gps_corr_abs",      &gps_corr_abs);
//    tree_final->SetBranchAddress("time_corr_abs",     &time_corr_abs);

    // Saving the original adress
//    uint32_t pps_cpt_corr, pps_cpt_init, ts_corr, gps_clock;
//    Double_t time_corr;
//    uint32_t pps_cpt_corr_noabs, pps_cpt_init_noabs, ts_corr_abs, gps_clock_abs;
//    Double_t time_corr_abs;

//    Finds the first entry where the absolute time gps is 97252 (number of second at 23rd june 23:00:52 with origin being 22nd june 20:00:00) return this entry
    uint32_t ref = 97252;
    uint32_t inc_index = tree_init->GetEntries();
    for (uint32_t i=0; i<tree_init->GetEntries(); i++) {
        tree_init->GetEntry(i);
        if (gps_corr - abs_gps_time_ref == ref){
            std::cout << "GPS absolute time of " << ref << " incremented at entry " << i << " - Returning this value." << std::endl;
            std::cout << "Absolute time chosen for aligning PPS and GPS : " << gps_corr - abs_gps_time_ref << std::endl;
            inc_index = i;
            break;
            }
    }
    if (inc_index == tree_init->GetEntries()){
        std::cerr << "ERROR : No entry found for the searched GPS absolute time." << std::endl;
    }

    tree_init->GetEntry(inc_index);
    uint32_t abs_pps_time_ref = pps_cpt_corr - gps_corr + abs_gps_time_ref;
    std::cout << "      Time correction applied to PSS  :  " << abs_pps_time_ref << "      Corresponding entry index : " << inc_index << std::endl;
    std::cout << "      Values of the PPS and GPS value  " << pps_cpt_corr - abs_pps_time_ref << "    " << gps_corr - abs_gps_time_ref << std::endl;
    uint64_t nentries = tree_init->GetEntries();
    for (uint64_t i = 0; i<nentries; i++) {
        tree_init->GetEntry(i);

        ts_init_noabs = ts_init;
        pps_cpt_init_noabs = pps_cpt_init;
        pps_cpt_corr_noabs = pps_cpt_corr;
        gps_corr_noabs = gps_corr;
        ts_corr_abs = ts_corr;
        pps_cpt_corr_abs = pps_cpt_corr - abs_pps_time_ref;
        pps_cpt_nocorr_abs = pps_cpt_init - abs_pps_time_ref;
        gps_corr_abs = gps_corr - abs_gps_time_ref;
        time_corr_abs = pps_cpt_corr_abs + ts_corr_abs/40.e6;

        tree_final->Fill();
    }
    file_init->Close();
    file_final->Write();
    file_final->Close();
}
