#include "CorrectionTimes.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TF1.h>
#include <Rtypes.h>
#include <string>
#include <fstream>
#include <array>
#include <TApplication.h>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include <TTreeIndex.h>
#include <TLeaf.h>

bool InWindow(const Double_t mid_interval, const Double_t half_window, const Double_t tested_value) {
    if (tested_value >= mid_interval - half_window && tested_value <= mid_interval + half_window) {
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

std::vector<size_t> sort_index(const std::vector<Double_t>& vec) {
    std::vector<size_t> index(vec.size());
    for (size_t i = 0; i < index.size(); i++) {
        index[i] = i;
    }
    std::sort(index.begin(), index.end(), [&vec](size_t i1, size_t i2) {
        return vec[i1] < vec[i2]; // To sort increasingly
    });
    return index;
}

////    !!! Voir these de Anne Meyer pour le fonctionnement des DSSSD

void ADCFormatingIJCLAB(const uint16_t *ener1, const uint16_t *ener2, Double_t &energy_adc, uint16_t *adc_channels, const std::string detector, const Double_t *correction, const uint16_t *pedestals) {
    //
    // Function used to transform the energy in ADC in different channels into keV energies
    // Size of ener_1 and ener_2 is not fixed as it might change with different detectors
    //
    // Method form D2B (maud):
    // Recoverall the ADC channels and substract the pedestal
    // Divide the value by the correction value (between 0 and 1)
    // Use the function to obtain energy : ax+b with a and b depending on temperature
    if (detector == "maud") {
        // Filling the adc_channel variable with pedestal-corrected ADC values
        for (uint64_t i = 0; i<32; i++){
            if (ener1[i] > pedestals[i]) {
                adc_channels[i] = ener1[i] - pedestals[i];
            } else {
                adc_channels[i] = 0;
            }
            if (ener2[i] > pedestals[i+32]) {
                adc_channels[i+32] = ener2[i] - pedestals[i+32];
            } else {
                adc_channels[i+32] = 0;
            }
        }
        // Searching the channel this highest ADC
        uint64_t ref_pixel_idx = 0;
        uint16_t max_adc = adc_channels[0];
        for (uint64_t i = 1; i<64; i++){
            if (adc_channels[i] > max_adc) {
                ref_pixel_idx = i;
                max_adc = adc_channels[i];
            }
        }
        // Calculation of the event's energy (in ADC) by summing the channels' ADC and applying a correction
        // The correction is due to the event being in the centre (corr ~ 1) or on the edge (corr < 1) of the crystal
        Double_t ener = 0;
        for (uint64_t i = 0; i<64; i++){
            ener += static_cast<Double_t>(adc_channels[i]);
        }
        energy_adc = ener / correction[ref_pixel_idx] / 64.;
    } else if (detector == "dssd") {
        // Filling the adc_channel variable ADC values
        energy_adc = 0.;
        for (uint64_t i = 0; i<32; i++){
            adc_channels[i] = ener1[i];
            adc_channels[i+32] = ener2[i];
        }
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud, dssd, cea or ucd" << std::endl;
    }
    return;
}

void ADCFormatingCEA(const Double_t *ener1, const Double_t *ener2, Double_t &energy_adc, uint16_t *adc_channels) {
    //
    // Function used to transform the energy in ADC in different channels into keV energies
    // Size of ener_1 and ener_2 is not fixed as it might change with different detectors
    //
    energy_adc = 0.;
    for (uint64_t i = 0; i<32; i++){
        adc_channels[i] = static_cast<uint16_t>(ener1[i]);
        adc_channels[i+32] = static_cast<uint16_t>(ener2[i]);
    }
    return;
}

Double_t ADCtoEnergyMaud(const Double_t ener_adc, const Double_t temp, const Double_t correction) {
    Double_t energy;
    // Calibration with lab measurements
//        Double_t alin = 0.0011086 * temp + 0.3957864;
//        Double_t blin = 0.157769 * temp - 16.251655;
//        energy = alin * energy_adc + blin;
    // Calibration with on flight data, made by Claire
//        Double_t slope = 0.0009 * temp + 0.4626;
//        energy = slope * ener_adc;
//        // Calibration with on flight data, made by Claire
//        energy =  1.4243E-5 * ener_adc * ener_adc + 0.386578 * ener_adc;
    // Calibration with correction
    energy =  correction * (1.4243E-5 * ener_adc * ener_adc + 0.386578 * ener_adc);
    return energy;
}

Double_t ADCtoEnergyDSSD(const uint16_t *adc_channels, Double_t *adc_channels_calib, Double_t *adc_channels_calib_corr, const Double_t temp, const Double_t coeffs[64][8]) {
    Double_t energy = 0.;
//    Double_t coef_a;
//    Double_t coef_b;
//    Double_t coef_c;
//    Double_t coef_d;
    Double_t coef_a_25;
    Double_t coef_b_25;
    Double_t coef_c_25;
    Double_t coef_d_25;
    Double_t temp_calib_file = 25.;
    std::vector<Double_t> a_corr = {-1.039, -0.855, -1.326, -1.08 , -0.989, -0.703, -0.523, -0.455, -0.955, -0.738, -0.312, -0.777, -0.556, -1.005, -0.86 , -0.581, -0.555, -0.506, -0.559, -0.605, -0.528, -0.635, -0.596, -0.308, -0.574, -0.218, -0.291, -0.494, -0.976, -0.391, -0.884, -0.727};
    std::vector<Double_t> b_corr = {554.288, 578.023, 670.033, 561.88 , 561.365, 558.516, 537.331, 563.017, 565.101, 542.656, 536.488, 548.783, 539.049, 542.394, 536.056, 548.699, 542.01 , 535.695, 541.572, 539.59 , 544.178, 524.454, 535.094, 531.449, 544.596, 537.483, 551.105, 546.186, 542.047, 555.329, 564.022, 570.375};
    Double_t corr;
    for (uint64_t i = 0; i<64; i++){
        // CAREFULL the first component of coeffs has as origin 25°C
        coef_a_25 = coeffs[i][0];
        coef_b_25 = coeffs[i][2];
        coef_c_25 = coeffs[i][4];
        coef_d_25 = coeffs[i][6];
        // We multiply by 1000 to have the correct unit (keV)
        adc_channels_calib[i] = 1000 * (coef_a_25 + coef_b_25 * adc_channels[i] + coef_c_25 * adc_channels[i] * adc_channels[i] + coef_d_25 * adc_channels[i] * adc_channels[i] * adc_channels[i]);
        if (adc_channels_calib[i] < 0.) {
            adc_channels_calib[i] = 0.;
        }

//        // CAREFULL the first component of coeffs has as origin 25°C
//        coef_a = coeffs[i][0] + coeffs[i][1] * (temp - temp_calib_file);
//        coef_b = coeffs[i][2] + coeffs[i][3] * (temp - temp_calib_file);
//        coef_c = coeffs[i][4] + coeffs[i][5] * (temp - temp_calib_file);
//        coef_d = coeffs[i][6] + coeffs[i][7] * (temp - temp_calib_file);
//        adc_channels_calib_corr[i] = 1000 * (coef_a + coef_b * adc_channels[i] + coef_c * adc_channels[i] * adc_channels[i] + coef_d * adc_channels[i] * adc_channels[i] * adc_channels[i]);
//        std::cout << adc_channels_calib_corr[i] << "  " << adc_channels_calib[i] << std::endl;
        // Latest correction, proportional correction with respect to 25° and to the bump
        // Bump variation is calculated, bump at 25° according to the variation is obtained
        // correction is obtained by comparing bump at T and at 25°
        // Then we have Ecorr = E / corr(T)
        if (i >= 32) {
//            corr = (a_corr[i - 32] * temp + b_corr[i - 32]) / (a_corr[i - 32] * temp_calib_file + b_corr[i - 32]);
            corr = (a_corr[i - 32] * temp + b_corr[i - 32]) / 676;
//            std::cout << << " " << << " " << corr << std::endl;
//            std::cout << corr << std::endl;
            adc_channels_calib_corr[i] = adc_channels_calib[i] / corr;
        } else {
            adc_channels_calib_corr[i] = adc_channels_calib[i];
        }
    }
    uint64_t max_idx = 0;
    for (uint64_t i = 1; i<64; i++){
        if (adc_channels_calib_corr[i] > adc_channels_calib_corr[i-1]) {
            max_idx = i;
        }
    }
    energy = adc_channels_calib_corr[max_idx];
    return energy;
}

Double_t ADCtoEnergyCEA(const uint16_t *adc_channels, Double_t *adc_channels_calib, const Double_t temp, const Double_t coeffs[64][8]) {
    Double_t energy = 0.;
    Double_t coef_a;
    Double_t coef_b;
    Double_t coef_c;
    Double_t coef_d;
    for (uint64_t i = 0; i<64; i++){
        coef_a = coeffs[i][0] + coeffs[i][1] * temp;
        coef_b = coeffs[i][2] + coeffs[i][3] * temp;
        coef_c = coeffs[i][4] + coeffs[i][5] * temp;
        coef_d = coeffs[i][6] + coeffs[i][7] * temp;
        adc_channels_calib[i] = 1000 * (coef_a + coef_b * adc_channels[i] + coef_c * adc_channels[i] * adc_channels[i] + coef_d * adc_channels[i] * adc_channels[i] * adc_channels[i]);

        if (adc_channels_calib[i] > 0.) {
            energy += adc_channels_calib[i];
        } else {
            adc_channels_calib[i] = 0.;
        }
    }
    return energy;
}

void ADCTempCorrectionDSSD(const uint16_t *adc_channels, Double_t *adc_channels_corr, const Double_t temp) {
    // a and b parameters of the affine function
    std::vector<Double_t> peak_a = {0.036716, 0.131445, 0.062599, 0.095734, 0.102531, 0.104413, 0.140915, 0.136512, 0.098898, 0.16447, 0.174816, 0.166083, 0.189503, 0.16458, 0.134001, 0.18001, 0.181169, 0.202537, 0.184703, 0.193654, 0.181026, 0.194065, 0.18256, 0.172182, 0.462903, 0.12925, 0.202558, 0.172908, 0.095815, 0.08862, 0.070518, 0.101594};
    std::vector<Double_t> peak_b = {90.271101, 76.100973, 76.508661, 141.380324, 74.719918, 73.90952, 137.309308, 72.493295, 74.396593, 73.718072, 72.328855, 71.554289, 71.762032, 73.482051, 72.323857, 71.98295, 72.441124, 71.740608, 72.07533, 71.566482, 73.003177, 71.484204, 71.522837, 71.590043, 130.925712, 238.847331, 73.02934, 72.792331, 145.007932, 258.649297, 201.89576, 84.206379};
    std::vector<Double_t> bulk_a = {-1.015557, -0.787783, -1.262267, -0.779096, -0.98369, -0.46709, -0.47222, -0.3871, -0.683492, -0.649886, -0.301839, -0.648335, -0.422, -0.627697, -0.632843, -0.542424, -0.494638, -0.489997, -0.354492, -0.711634, -0.504241, -0.661699, -0.501675, -0.374338, -0.351215, -0.183012, -0.211629, -0.498082, -0.676097, -0.409012, -1.041214, -0.5786};
    std::vector<Double_t> bulk_b = {575.174915, 610.581457, 590.756867, 573.391298, 591.800576, 579.43737, 569.640771, 577.149503, 587.798512, 584.149709, 573.035371, 569.42912, 569.633759, 584.4602, 565.89034, 577.212994, 580.591879, 570.205355, 572.179833, 574.642222, 576.478365, 567.970641, 567.467592, 564.339591, 561.814468, 557.396933, 582.679154, 574.200107, 570.794763, 567.394666, 584.253568, 581.42604};

    std::vector<std::vector<Double_t>> peak_adc = {
        {90.6881, 78.1919, 77.8653, 143.262, 76.5146, 75.6742, 140.12, 74.6548, 76.2541, 76.4875, 75.1697, 74.4551, 74.7664, 76.3069, 74.7252, 74.9019, 75.4198, 74.94, 75.1563, 74.6333, 76.0561, 74.5017, 74.3554, 74.2974, 137.674, 247.856, 76.2419, 75.631, 147.43, 264.252, 203.583, 85.1128},
        {90.57, 77.1416, 77.668, 141.743, 75.7145, 74.991, 138.108, 73.698, 75.2459, 74.8408, 73.4489, 72.4524, 72.8914, 74.5194, 73.2314, 73.0635, 73.4874, 72.9512, 73.1201, 72.7235, 74.0617, 72.6125, 72.6595, 72.5875, 133.888, 234.842, 74.2762, 73.8409, 144.966, 255.983, 201.791, 85.6977},
        {90.9995, 76.1788, 75.962, 141.226, 74.6853, 73.4866, 137.142, 72.5588, 74.3386, 73.7731, 72.4255, 71.7156, 71.8925, 73.5476, 72.356, 72.189, 72.5503, 71.9745, 72.0912, 71.7015, 73.1307, 71.6362, 71.672, 71.9685, 131.468, 231.438, 73.2042, 72.8193, 144.455, 246.698, 201.248, 85.2755},
        {89.7645, 75.4348, 75.448, 141.042, 74.0429, 72.7642, 136.739, 71.4654, 73.7001, 72.8441, 71.4768, 70.7566, 70.7234, 72.6347, 71.6124, 71.1803, 71.6746, 70.8928, 71.2551, 70.7174, 72.1329, 70.6096, 70.864, 70.4891, 130.188, 239.726, 71.995, 71.8133, 144.38, 257.287, 201.445, 83.9695},
        {89.7303, 74.8247, 74.5934, 140.731, 73.516, 71.5873, 136.298, 70.7879, 73.2497, 72.1315, 70.7434, 69.9521, 70.2194, 72.0507, 71.0668, 70.487, 70.9175, 70.0817, 70.5895, 70.0644, 71.4679, 70.0643, 69.8823, 69.8783, 128.24, 241.254, 71.2397, 71.1628, 144.368, 266.332, 201.429, 83.1174},
        {89.8969, 74.6727, 76.6052, 140.26, 73.6671, 74.3619, 135.797, 71.505, 73.3877, 71.981, 70.5166, 69.8374, 69.9852, 71.6267, 70.7812, 69.9119, 70.4011, 69.4496, 70.0035, 69.39, 70.9461, 69.4047, 69.4954, 69.8965, 124.349, 237.261, 70.9427, 71.0418, 144.145, 259.125, 201.261, 82.7058}
    };

    std::vector<std::vector<Double_t>> bulk_adc = {
        {557.566, 597.997, 570.104, 561.902, 572.592, 571.81, 560.698, 564.902, 577.489, 575.975, 563.334, 551.994, 562.456, 577.246, 560.419, 568.872, 573.038, 561.024, 571.325, 564.126, 568.062, 558.912, 552.399, 557.039, 549.489, 551.144, 576.482, 563.05, 561.278, 558.804, 567.102, 571.466},
        {567.124, 607.184, 584.25, 564.591, 586.862, 576.591, 568.159, 578.319, 581.216, 575.968, 570.129, 568.391, 564.915, 578.611, 560.316, 573.044, 577.639, 568.38, 567.366, 567.335, 575.128, 566.384, 564.694, 561.687, 562.314, 556.455, 582.583, 571.691, 564.398, 563.34, 577.824, 578.473},
        {577.552, 605.255, 582.335, 575.562, 588.193, 580.338, 565.023, 574.373, 592.037, 584.452, 575.433, 565.617, 565.535, 583.684, 559.379, 572.968, 577.012, 565.329, 566.038, 577.297, 574.193, 559.948, 568.789, 564.763, 555.852, 558.849, 579.779, 571.832, 568.674, 571.789, 579.991, 576.884},
        {575.727, 610.195, 595.82, 576.242, 598.649, 580.954, 573.912, 579.731, 587.524, 588.091, 575.836, 573.634, 577.043, 585.668, 564.369, 584.861, 586.256, 573.501, 575.1, 576.88, 573.712, 566.595, 573.038, 566.131, 564.437, 560.755, 586.95, 577.668, 569.882, 565.453, 592.644, 587.238},
        {583.992, 622.168, 594.827, 580.657, 603.231, 581.755, 571.914, 583.017, 594.566, 588.984, 577.451, 572.019, 573.215, 586.812, 572.397, 579.549, 581.471, 575.038, 572.039, 580.899, 578.067, 572.33, 569.481, 567.506, 568.859, 559.371, 586.73, 582.554, 580.716, 571.619, 584.588, 587.59},
        {587.595, 618.91, 608.909, 582.329, 602.151, 585.741, 575.807, 580.956, 596.215, 592.131, 575.036, 578.206, 573.377, 593.885, 574.711, 583.793, 587.282, 575.743, 578.279, 583.472, 584.176, 577.523, 573.344, 568.363, 563.685, 558.673, 583.518, 577.16, 577.321, 572.408, 598.239, 586.734}
    };
    size_t temperature_ite;

    if (temp >= 10) {
        temperature_ite = 0;
    } else if (temp < 10 && temp >= 5) {
        temperature_ite = 1;
    } else if (temp < 5 && temp >= 0) {
        temperature_ite = 2;
    } else if (temp < 0 && temp >= -5) {
        temperature_ite = 3;
    } else if (temp < -5 && temp >= -10) {
        temperature_ite = 4;
    } else if (temp < -10) {
        temperature_ite = 5;
    }
    Double_t peak_adc_temp;
    Double_t bulk_adc_temp;
    Double_t peak_adc_25;
    Double_t bulk_adc_25;
    Double_t corr_peak;
    Double_t corr_bulk;
    Double_t corr_a;
    Double_t corr_b;
//    std::cout << "temp : " << temp << std::endl;
    for (size_t i = 0; i < 32; i++) {
        adc_channels_corr[i] = adc_channels[i];
        // Correcting the n face channels if it's != 0
        if (adc_channels[i + 32] == 0) {
            adc_channels_corr[i + 32] = adc_channels[i + 32];
        } else {
//            std::cout << "chan : " << i << std::endl;
            // Peak and bulk ADC at T temp and 25°, obtained with affine function from fitting raw spectra
            peak_adc_temp = peak_b[i] + peak_a[i] * temp;
            bulk_adc_temp = bulk_b[i] + bulk_a[i] * temp;
//            peak_adc_temp = peak_adc[temperature_ite][i];
//            bulk_adc_temp = bulk_adc[temperature_ite][i];
            peak_adc_25 = peak_b[i] + peak_a[i] * 25;
            bulk_adc_25 = bulk_b[i] + bulk_a[i] * 25;
//            std::cout << "peak_adc_temp : " << peak_adc_temp << std::endl;
//            std::cout << "peak_adc_25 : " << peak_adc_25 << std::endl;
//            std::cout << "bulk_adc_temp : " << bulk_adc_temp << std::endl;
//            std::cout << "bulk_adc_25 : " << bulk_adc_25 << std::endl;
//            // Correction applicable to peak and bulk
//            corr_peak = peak_adc_25 / peak_adc_temp;
//            corr_bulk = bulk_adc_25 / bulk_adc_temp;
////            std::cout << "corr_peak : " << corr_peak << std::endl;
////            std::cout << "corr_bulk : " << corr_bulk << std::endl;
            // Using these 2 correction points to estimate an affine correction function corr = f(ADC)
//            corr_a = (corr_bulk - corr_peak) / (bulk_adc_temp - peak_adc_temp);
//            corr_b = corr_bulk - corr_a * bulk_adc_temp;
            corr_a = (bulk_adc_25 - peak_adc_25) / (bulk_adc_temp - peak_adc_temp);
            corr_b = bulk_adc_25 - corr_a * bulk_adc_temp;
//            std::cout << "corr slope : " << corr_a << std::endl;
//            std::cout << "corr b origin : " << corr_b << std::endl;
//            adc_channels_corr[i + 32] = adc_channels[i + 32] * (corr_a * adc_channels[i + 32] + corr_b);
            adc_channels_corr[i + 32] = corr_a * adc_channels[i + 32] + corr_b;
//            std::cout << "adc_channels : " << adc_channels[i + 32] << std::endl;
//            std::cout << "adc_channels_corr : " << adc_channels_corr[i + 32] << std::endl;
        }
    }
//    std::exit(EXIT_SUCCESS);
}

void ADCTempCorrectionCEA(const uint16_t *adc_channels, Double_t *adc_channels_corr, const Double_t temp) {
    std::vector<Double_t> peak_a = {0.036716, 0.131445, 0.062599, 0.095734, 0.102531, 0.104413, 0.140915, 0.136512, 0.098898, 0.16447, 0.174816, 0.166083, 0.189503, 0.16458, 0.134001, 0.18001, 0.181169, 0.202537, 0.184703, 0.193654, 0.181026, 0.194065, 0.18256, 0.172182, 0.462903, 0.12925, 0.202558, 0.172908, 0.095815, 0.08862, 0.070518, 0.101594};
    std::vector<Double_t> peak_b = {90.271101, 76.100973, 76.508661, 141.380324, 74.719918, 73.90952, 137.309308, 72.493295, 74.396593, 73.718072, 72.328855, 71.554289, 71.762032, 73.482051, 72.323857, 71.98295, 72.441124, 71.740608, 72.07533, 71.566482, 73.003177, 71.484204, 71.522837, 71.590043, 130.925712, 238.847331, 73.02934, 72.792331, 145.007932, 258.649297, 201.89576, 84.206379};
    std::vector<Double_t> bulk_a = {-1.015557, -0.787783, -1.262267, -0.779096, -0.98369, -0.46709, -0.47222, -0.3871, -0.683492, -0.649886, -0.301839, -0.648335, -0.422, -0.627697, -0.632843, -0.542424, -0.494638, -0.489997, -0.354492, -0.711634, -0.504241, -0.661699, -0.501675, -0.374338, -0.351215, -0.183012, -0.211629, -0.498082, -0.676097, -0.409012, -1.041214, -0.5786};
    std::vector<Double_t> bulk_b = {575.174915, 610.581457, 590.756867, 573.391298, 591.800576, 579.43737, 569.640771, 577.149503, 587.798512, 584.149709, 573.035371, 569.42912, 569.633759, 584.4602, 565.89034, 577.212994, 580.591879, 570.205355, 572.179833, 574.642222, 576.478365, 567.970641, 567.467592, 564.339591, 561.814468, 557.396933, 582.679154, 574.200107, 570.794763, 567.394666, 584.253568, 581.42604};

    std::vector<std::vector<Double_t>> peak_adc = {
        {90.6881, 78.1919, 77.8653, 143.262, 76.5146, 75.6742, 140.12, 74.6548, 76.2541, 76.4875, 75.1697, 74.4551, 74.7664, 76.3069, 74.7252, 74.9019, 75.4198, 74.94, 75.1563, 74.6333, 76.0561, 74.5017, 74.3554, 74.2974, 137.674, 247.856, 76.2419, 75.631, 147.43, 264.252, 203.583, 85.1128},
        {90.57, 77.1416, 77.668, 141.743, 75.7145, 74.991, 138.108, 73.698, 75.2459, 74.8408, 73.4489, 72.4524, 72.8914, 74.5194, 73.2314, 73.0635, 73.4874, 72.9512, 73.1201, 72.7235, 74.0617, 72.6125, 72.6595, 72.5875, 133.888, 234.842, 74.2762, 73.8409, 144.966, 255.983, 201.791, 85.6977},
        {90.9995, 76.1788, 75.962, 141.226, 74.6853, 73.4866, 137.142, 72.5588, 74.3386, 73.7731, 72.4255, 71.7156, 71.8925, 73.5476, 72.356, 72.189, 72.5503, 71.9745, 72.0912, 71.7015, 73.1307, 71.6362, 71.672, 71.9685, 131.468, 231.438, 73.2042, 72.8193, 144.455, 246.698, 201.248, 85.2755},
        {89.7645, 75.4348, 75.448, 141.042, 74.0429, 72.7642, 136.739, 71.4654, 73.7001, 72.8441, 71.4768, 70.7566, 70.7234, 72.6347, 71.6124, 71.1803, 71.6746, 70.8928, 71.2551, 70.7174, 72.1329, 70.6096, 70.864, 70.4891, 130.188, 239.726, 71.995, 71.8133, 144.38, 257.287, 201.445, 83.9695},
        {89.7303, 74.8247, 74.5934, 140.731, 73.516, 71.5873, 136.298, 70.7879, 73.2497, 72.1315, 70.7434, 69.9521, 70.2194, 72.0507, 71.0668, 70.487, 70.9175, 70.0817, 70.5895, 70.0644, 71.4679, 70.0643, 69.8823, 69.8783, 128.24, 241.254, 71.2397, 71.1628, 144.368, 266.332, 201.429, 83.1174},
        {89.8969, 74.6727, 76.6052, 140.26, 73.6671, 74.3619, 135.797, 71.505, 73.3877, 71.981, 70.5166, 69.8374, 69.9852, 71.6267, 70.7812, 69.9119, 70.4011, 69.4496, 70.0035, 69.39, 70.9461, 69.4047, 69.4954, 69.8965, 124.349, 237.261, 70.9427, 71.0418, 144.145, 259.125, 201.261, 82.7058}
    };

    std::vector<std::vector<Double_t>> bulk_adc = {
        {557.566, 597.997, 570.104, 561.902, 572.592, 571.81, 560.698, 564.902, 577.489, 575.975, 563.334, 551.994, 562.456, 577.246, 560.419, 568.872, 573.038, 561.024, 571.325, 564.126, 568.062, 558.912, 552.399, 557.039, 549.489, 551.144, 576.482, 563.05, 561.278, 558.804, 567.102, 571.466},
        {567.124, 607.184, 584.25, 564.591, 586.862, 576.591, 568.159, 578.319, 581.216, 575.968, 570.129, 568.391, 564.915, 578.611, 560.316, 573.044, 577.639, 568.38, 567.366, 567.335, 575.128, 566.384, 564.694, 561.687, 562.314, 556.455, 582.583, 571.691, 564.398, 563.34, 577.824, 578.473},
        {577.552, 605.255, 582.335, 575.562, 588.193, 580.338, 565.023, 574.373, 592.037, 584.452, 575.433, 565.617, 565.535, 583.684, 559.379, 572.968, 577.012, 565.329, 566.038, 577.297, 574.193, 559.948, 568.789, 564.763, 555.852, 558.849, 579.779, 571.832, 568.674, 571.789, 579.991, 576.884},
        {575.727, 610.195, 595.82, 576.242, 598.649, 580.954, 573.912, 579.731, 587.524, 588.091, 575.836, 573.634, 577.043, 585.668, 564.369, 584.861, 586.256, 573.501, 575.1, 576.88, 573.712, 566.595, 573.038, 566.131, 564.437, 560.755, 586.95, 577.668, 569.882, 565.453, 592.644, 587.238},
        {583.992, 622.168, 594.827, 580.657, 603.231, 581.755, 571.914, 583.017, 594.566, 588.984, 577.451, 572.019, 573.215, 586.812, 572.397, 579.549, 581.471, 575.038, 572.039, 580.899, 578.067, 572.33, 569.481, 567.506, 568.859, 559.371, 586.73, 582.554, 580.716, 571.619, 584.588, 587.59},
        {587.595, 618.91, 608.909, 582.329, 602.151, 585.741, 575.807, 580.956, 596.215, 592.131, 575.036, 578.206, 573.377, 593.885, 574.711, 583.793, 587.282, 575.743, 578.279, 583.472, 584.176, 577.523, 573.344, 568.363, 563.685, 558.673, 583.518, 577.16, 577.321, 572.408, 598.239, 586.734}
    };
    size_t temperature_ite;

    if (temp >= 10) {
        temperature_ite = 0;
    } else if (temp < 10 && temp >= 5) {
        temperature_ite = 1;
    } else if (temp < 5 && temp >= 0) {
        temperature_ite = 2;
    } else if (temp < 0 && temp >= -5) {
        temperature_ite = 3;
    } else if (temp < -5 && temp >= -10) {
        temperature_ite = 4;
    } else if (temp < -10) {
        temperature_ite = 5;
    }
    Double_t peak_adc_temp;
    Double_t bulk_adc_temp;
    Double_t peak_adc_25;
    Double_t bulk_adc_25;
    Double_t corr_peak;
    Double_t corr_bulk;
    Double_t corr_a;
    Double_t corr_b;
//    std::cout << "temp : " << temp << std::endl;
    for (size_t i = 0; i < 32; i++) {
        adc_channels_corr[i] = adc_channels[i];
        // Correcting the n face channels if it's != 0
        if (adc_channels[i + 32] == 0) {
            adc_channels_corr[i + 32] = adc_channels[i + 32];
        } else {
//            std::cout << "chan : " << i << std::endl;
            // Peak and bulk ADC at T temp and 25°, obtained with affine function from fitting raw spectra
            peak_adc_temp = peak_b[i] + peak_a[i] * temp;
            bulk_adc_temp = bulk_b[i] + bulk_a[i] * temp;
//            peak_adc_temp = peak_adc[temperature_ite][i];
//            bulk_adc_temp = bulk_adc[temperature_ite][i];
            peak_adc_25 = peak_b[i] + peak_a[i] * 25;
            bulk_adc_25 = bulk_b[i] + bulk_a[i] * 25;
//            std::cout << "peak_adc_temp : " << peak_adc_temp << std::endl;
//            std::cout << "peak_adc_25 : " << peak_adc_25 << std::endl;
//            std::cout << "bulk_adc_temp : " << bulk_adc_temp << std::endl;
//            std::cout << "bulk_adc_25 : " << bulk_adc_25 << std::endl;
//            // Correction applicable to peak and bulk
//            corr_peak = peak_adc_25 / peak_adc_temp;
//            corr_bulk = bulk_adc_25 / bulk_adc_temp;
////            std::cout << "corr_peak : " << corr_peak << std::endl;
////            std::cout << "corr_bulk : " << corr_bulk << std::endl;
            // Using these 2 correction points to estimate an affine correction function corr = f(ADC)
//            corr_a = (corr_bulk - corr_peak) / (bulk_adc_temp - peak_adc_temp);
//            corr_b = corr_bulk - corr_a * bulk_adc_temp;
            corr_a = (bulk_adc_25 - peak_adc_25) / (bulk_adc_temp - peak_adc_temp);
            corr_b = bulk_adc_25 - corr_a * bulk_adc_temp;
//            std::cout << "corr slope : " << corr_a << std::endl;
//            std::cout << "corr b origin : " << corr_b << std::endl;
//            adc_channels_corr[i + 32] = adc_channels[i + 32] * (corr_a * adc_channels[i + 32] + corr_b);
            adc_channels_corr[i + 32] = corr_a * adc_channels[i + 32] + corr_b;
//            std::cout << "adc_channels : " << adc_channels[i + 32] << std::endl;
//            std::cout << "adc_channels_corr : " << adc_channels_corr[i + 32] << std::endl;
        }
    }
//    std::exit(EXIT_SUCCESS);
}

void ExtractTemperatureValues(std::vector<Double_t> &tempvec, std::vector<int32_t> &tsvec, const std::string timestamp_filename) {
    size_t start, end;
    uint16_t tokencounting;
    std::ifstream timestamp_file(timestamp_filename); // Open file
    if (!timestamp_file.is_open()) {
        std::cerr << "Error : Opening of the file failed !" << std::endl;
        return; //
    }
    std::string tsline;
    const std::string delim = "\t";
    if (std::getline(timestamp_file, tsline)) {
        std::cout << "Reading time vs temperature file" << std::endl;
        std::cout << "First line giving columns : " << tsline << std::endl;
    } else {
        std::cerr << "Unable to get the first line of the time vs temperature file" << std::endl;
    }
    while (std::getline(timestamp_file, tsline)) {
        tokencounting = 0;
        start = 0;
        end = 0;

        while ((end = tsline.find_first_of(delim, start)) != std::string::npos) {
            if (end != start) { // Ignore consecutive delimitors
                if (tokencounting == 2) {
                    tsvec.push_back(static_cast<int32_t>(std::stoi(tsline.substr(start, end - start))));
                }
                tokencounting += 1;
            }
            start = end + 1;
        }
        if (start < tsline.size()) {
            tempvec.push_back(std::stod(tsline.substr(start)));
        }
    }
    timestamp_file.close();
}

void FindTemperature(const int32_t pps_timestamp, Double_t &temp_final, uint64_t &ite_init, const std::vector<int32_t> timestamp_vec, std::vector<Double_t> temperature_vec) {
    uint64_t loop_starter = ite_init;
    for (uint64_t i = loop_starter; i<timestamp_vec.size(); i++) {
        if (pps_timestamp <= timestamp_vec[i]) {
            if (pps_timestamp == timestamp_vec[i]) {
                temp_final = temperature_vec[i];
            } else {
                if (i == 0) {
                    temp_final = temperature_vec[i];
                } else {
                    temp_final = temperature_vec[i - 1] + (temperature_vec[i] - temperature_vec[i - 1]) / (timestamp_vec[i] - timestamp_vec[i - 1]) * (pps_timestamp - timestamp_vec[i - 1]);
                }
            }
            return;
        }
        ite_init = i;
    }
    return;
}

void ExtractCorrectionValues(std::vector<Double_t> &corrvec, std::vector<int32_t> &tsvec, const std::string correction_filename) {
    size_t start, end;
    uint16_t tokencounting;
    std::ifstream timestamp_file(correction_filename); // Open file
    if (!timestamp_file.is_open()) {
        std::cerr << "Error : Opening of the file failed !" << std::endl;
        return; //
    }
    std::string tsline;
    const std::string delim = " ";
    if (std::getline(timestamp_file, tsline)) {
        std::cout << "Reading time vs correction file" << std::endl;
        std::cout << "First line giving columns : " << tsline << std::endl;
    } else {
        std::cerr << "Unable to get the first line of the time vs temperature file" << std::endl;
    }
    while (std::getline(timestamp_file, tsline)) {
        tokencounting = 0;
        start = 0;
        end = 0;

        while ((end = tsline.find_first_of(delim, start)) != std::string::npos) {
            if (end != start) { // Ignore consecutive delimitors
                if (tokencounting == 0) {
                    tsvec.push_back(static_cast<int32_t>(std::stoi(tsline.substr(start, end - start))));
                }
                tokencounting += 1;
            }
            start = end + 1;
        }
        if (start < tsline.size()) {
            corrvec.push_back(std::stod(tsline.substr(start)));
        }
    }
    timestamp_file.close();
}

void FindCorrection(const int32_t pps_timestamp, Double_t &correction, uint64_t &ite_init, const std::vector<int32_t> timestamp_vec, std::vector<Double_t> correction_vec) {
    uint64_t loop_starter = ite_init;
    for (uint64_t i = loop_starter; i<timestamp_vec.size(); i++) {
        if (pps_timestamp == timestamp_vec[i]) {
            correction = correction_vec[i];
            return;
        }
        ite_init = i;
    }
    return;
}

void ExtractDSSDCoefs(Double_t coeffs[64][8], const std::string dssd_p_coefs_filename, const std::string dssd_n_coefs_filename) {
    size_t start, end;
    uint16_t tokencounting;
    size_t file_idx = 0;

    std::ifstream p_coef_file(dssd_p_coefs_filename); // Open p face file
    if (!p_coef_file.is_open()) {
        std::cerr << "Error : Opening of the file failed !" << std::endl;
        return; //
    }

    std::ifstream n_coef_file(dssd_n_coefs_filename); // Open p face file
    if (!n_coef_file.is_open()) {
        std::cerr << "Error : Opening of the file failed !" << std::endl;
        return; //
    }

    std::string line;
    const std::string delim = "\t";
//    a tester
    while (std::getline(p_coef_file, line)) {
        tokencounting = 0;
        start = 0;
        end = 0;
        while ((end = line.find_first_of(delim, start)) != std::string::npos) {
            if (end != start) { // Ignore consecutive delimitors
                if (tokencounting == 1) {
                    coeffs[file_idx][0] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 2) {
                    coeffs[file_idx][1] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 3) {
                    coeffs[file_idx][2] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 4) {
                    coeffs[file_idx][3] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 5) {
                    coeffs[file_idx][4] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 6) {
                    coeffs[file_idx][5] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 7) {
                    coeffs[file_idx][6] = std::stod(line.substr(start, end - start));
                }
                tokencounting += 1;
            }
            start = end + 1;
        }
        if (start < line.size()) {
            coeffs[file_idx][7] = std::stod(line.substr(start, end - start));
        }
        file_idx += 1;
    }
    if (file_idx != 32) {
        std::cerr << "ERROR : The index is supposed to be equal to 32. Got " << file_idx << " instead." << std::endl;
    }
    while (std::getline(n_coef_file, line)) {
        tokencounting = 0;
        start = 0;
        end = 0;
        while ((end = line.find_first_of(delim, start)) != std::string::npos) {
            if (end != start) { // Ignore consecutive delimitors
                if (tokencounting == 1) {
                    coeffs[file_idx][0] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 2) {
                    coeffs[file_idx][1] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 3) {
                    coeffs[file_idx][2] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 4) {
                    coeffs[file_idx][3] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 5) {
                    coeffs[file_idx][4] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 6) {
                    coeffs[file_idx][5] = std::stod(line.substr(start, end - start));
                } else if (tokencounting == 7) {
                    coeffs[file_idx][6] = std::stod(line.substr(start, end - start));
                }
                tokencounting += 1;
            }
            start = end + 1;
        }
        if (start < line.size()) {
            coeffs[file_idx][7] = std::stod(line.substr(start, end - start));
        }
        file_idx += 1;
    }
    p_coef_file.close();
    n_coef_file.close();
}

void FindPosition(uint16_t *adc_channels, Double_t *position, std::string detector) {
    if (detector == "maud") {
        static const std::map<int32_t, std::pair<int32_t, int32_t>> geometry {
        { 0, {5,3}}, { 1, {6,3}}, { 2, {7,3}}, { 3, {7,2}}, { 4, {7,1}}, { 5, {6,2}}, { 6, {6,1}}, { 7, {7,0}},
        { 8, {5,0}}, { 9, {6,0}}, {10, {4,0}}, {11, {5,2}}, {12, {5,1}}, {13, {4,3}}, {14, {4,2}}, {15, {4,1}},
        {16, {4,7}}, {17, {4,6}}, {18, {5,7}}, {19, {4,5}}, {20, {4,4}}, {21, {5,6}}, {22, {7,7}}, {23, {6,7}},
        {24, {7,6}}, {25, {5,5}}, {26, {7,5}}, {27, {6,6}}, {28, {7,4}}, {29, {6,5}}, {30, {6,4}}, {31, {5,4}},
        {32, {2,4}}, {33, {1,5}}, {34, {1,4}}, {35, {0,4}}, {36, {1,6}}, {37, {0,5}}, {38, {2,6}}, {39, {0,6}},
        {40, {2,5}}, {41, {0,7}}, {42, {2,7}}, {43, {1,7}}, {44, {3,4}}, {45, {3,6}}, {46, {3,5}}, {47, {3,7}},
        {48, {3,1}}, {49, {2,0}}, {50, {1,3}}, {51, {1,2}}, {52, {1,0}}, {53, {2,3}}, {54, {0,0}}, {55, {3,2}},
        {56, {0,1}}, {57, {2,2}}, {58, {0,2}}, {59, {3,3}}, {60, {0,3}}, {61, {3,0}}, {62, {1,1}}, {63, {2,1}}
        };

        // CARE the disposition remains pretty hard to precise : with this disposition y is the thickness, x goes from left to right and z from top to bottom
        static const size_t PixIDx[8][8] {
        {54, 56, 58, 60, 35, 37, 39, 41},
        {52, 62, 51, 50, 34, 33, 36, 43},
        {49, 63, 57, 53, 32, 40, 38, 42},
        {61, 48, 55, 59, 44, 46, 45, 47},
        {10, 15, 14, 13, 20, 19, 17, 16},
        { 8, 12, 11,  0, 31, 25, 21, 18},
        { 9,  6,  5,  1, 30, 29, 27, 23},
        { 7,  4,  3,  2, 28, 26, 24, 22}
        };

        uint64_t max_idx = 0;
        for (uint64_t i = 1; i<64; i++){
            if (adc_channels[i] > adc_channels[i-1]) {
                max_idx = i;
            }
        }
        // CARE MAKE SUR THE
        // Position are given with respect to the centre of the detector
        // Centroid method is used
        position[0] = 0.;
        position[1] = 0.;
        position[2] = 0.;
        uint32_t adcsum = 0;
//        std::cout << "A tester !!!" << std::endl;
        for (int64_t iidx : {geometry.at(max_idx).first - 1, geometry.at(max_idx).first + 1}) {
            for (int64_t jidx : {geometry.at(max_idx).second - 1, geometry.at(max_idx).second + 1}) {
                if (iidx >= 0 && iidx < 8 && jidx >= 0 && jidx < 8) {
                    position[0] += (jidx * 0.6375 - 2.23125) * static_cast<Double_t>(adc_channels[PixIDx[iidx][jidx]]);
                    position[2] += (-iidx * 0.6375 + 2.23125) * static_cast<Double_t>(adc_channels[PixIDx[iidx][jidx]]);
                    adcsum += static_cast<uint16_t>(adc_channels[PixIDx[iidx][jidx]]);
                }
            }
        }
        // Applying the offset (position of the detector in the instrument) and dividing by adcsum for the ceintroid
        Double_t xoffset = 0.;
        Double_t yoffset = 0.;
        Double_t zoffset = 0.;
        position[0] = xoffset + position[0] / static_cast<Double_t>(adcsum);
        position[1] = yoffset + 0.5;
        position[2] = zoffset + position[2] / static_cast<Double_t>(adcsum);
    } else if (detector == "dssd") {
        position[0] = 0.;
        position[1] = 0.;
        position[2] = 0.;
    } else if (detector == "cea") {
        position[0] = 0.;
        position[1] = 0.;
        position[2] = 0.;
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud, dssd, cea or ucd" << std::endl;
    }
    return;
}

void GetUCDPosition(Double_t position_final[3], const Short_t position_init) {
    static const std::map<int32_t, std::array<int32_t, 3>> geometry {
    { 0, {0, 1, 2}}, { 1, {0, 1, 2}}, { 2, {0, 1, 2}}, { 3, {0, 1, 2}},
    { 4, {0, 1, 2}}, { 5, {0, 1, 2}}, { 6, {0, 1, 2}}, { 7, {0, 1, 2}},
    { 8, {0, 1, 2}}, { 9, {0, 1, 2}}, {10, {0, 1, 2}}, {11, {0, 1, 2}},
    {12, {0, 1, 2}}, {13, {0, 1, 2}}, {14, {0, 1, 2}}, {15, {0, 1, 2}},
    };

    Double_t xoffset = 0.;
    Double_t yoffset = 0.;
    Double_t zoffset = 0.;

    position_final[0] = geometry.at(position_init - 1)[0] + xoffset;
    position_final[1] = geometry.at(position_init - 1)[1] + yoffset;
    position_final[2] = geometry.at(position_init - 1)[2] + zoffset;
}

void CombineDSSSDEvents(const std::string &init_tree_name, const std::string &final_tree_name) {
    //
    // Function used to combine the DSSSD channels for the 2 faces (P and N)
    //
    // Declaring file_final
    TFile* file_final = nullptr;
    // Declaring tree_final
    TTree* tree_final = nullptr;

    // Extraction of the root tree
    TChain *tree_init = new TChain("Events");
    tree_init->Add(init_tree_name.c_str());

    uint32_t timestamp_init;
    uint32_t pps_cpt_init;
    uint32_t pps_info_init;
    uint32_t event_cpt_init;
    uint8_t chain_init;
    uint8_t coincidence_init;
    uint8_t chip_data_init;
    uint8_t analog_trigger_init;
    uint8_t seu_init;
    uint32_t ch_status_init;
    uint16_t ref_channel_init;
    uint16_t sample_init[32];
    uint16_t cm_data_init;
    tree_init->SetBranchAddress("timestamp", &timestamp_init);
    tree_init->SetBranchAddress("pps_cpt", &pps_cpt_init);
    tree_init->SetBranchAddress("pps_info", &pps_info_init);
    tree_init->SetBranchAddress("event_cpt", &event_cpt_init);
    tree_init->SetBranchAddress("chain", &chain_init);
    tree_init->SetBranchAddress("coincidence", &coincidence_init);
    tree_init->SetBranchAddress("chip_data", &chip_data_init);
    tree_init->SetBranchAddress("analog_trigger", &analog_trigger_init);
    tree_init->SetBranchAddress("seu", &seu_init);
    tree_init->SetBranchAddress("ch_status", &ch_status_init);
    tree_init->SetBranchAddress("ref_channel", &ref_channel_init);
    tree_init->SetBranchAddress("sample", &sample_init);
    tree_init->SetBranchAddress("cm_data", &cm_data_init);

    file_final = new TFile(final_tree_name.c_str(), "RECREATE", "DSSSD combined file");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
    tree_final = new TTree("Events", "DSSSD combined file");

    uint32_t timestamp_final;
    uint32_t pps_cpt_final;
    uint32_t pps_info_final;
    uint32_t event_cpt_final;
//    uint8_t chain_final;
    uint8_t coincidence_final;
    uint8_t chip_data_final;
    uint8_t analog_trigger_final;
    uint8_t seu_final;
    uint32_t ch_status_final;
    uint16_t ref_channel_final;
    uint16_t sample_p[32];
    uint16_t sample_n[32];
    uint16_t cm_data_final;

    tree_final->Branch("timestamp", &timestamp_final, "timestamp/i");  // i = uint32_t
    tree_final->Branch("pps_cpt", &pps_cpt_final, "pps_cpt/i");        // i = uint32_t
    tree_final->Branch("pps_info", &pps_info_final, "pps_info/i");     // i = uint32_t
    tree_final->Branch("event_cpt", &event_cpt_final, "event_cpt/i");  // i = uint32_t
//    tree_final->Branch("chain", &chain_final, "chain/b");              // b = uint8_t
    tree_final->Branch("coincidence", &coincidence_final, "coincidence/b");
    tree_final->Branch("chip_data", &chip_data_final, "chip_data/b");
    tree_final->Branch("analog_trigger", &analog_trigger_final, "analog_trigger/b");
    tree_final->Branch("seu", &seu_final, "seu/b");
    tree_final->Branch("ch_status", &ch_status_final, "ch_status/i");
    tree_final->Branch("ref_channel", &ref_channel_final, "ref_channel/s");
    tree_final->Branch("sample_p", &sample_p, "sample_p[32]/s");
    tree_final->Branch("sample_n", &sample_n, "sample_n[32]/s");
    tree_final->Branch("cm_data", &cm_data_final, "cm_data/s");

    uint64_t nentries = tree_init->GetEntries();
    uint64_t ite_entries = 0;
    uint32_t temp_ts;
    uint32_t temp_pps;
    uint8_t temp_chain;
    while (ite_entries < nentries - 1) {
        // Getting ts, pps and chain for iteration ite_entries
        tree_init->GetEntry(ite_entries);
        temp_ts = timestamp_init;
        temp_pps = pps_cpt_init;
        temp_chain = chain_init;
        // Getting ts, pps and chain for iteration ite_entries + 1
        tree_init->GetEntry(ite_entries + 1);
        // Checking if the 2 events are the same but for p and n faces
        if (temp_ts == timestamp_init && temp_pps == pps_cpt_init) {
            if (temp_chain != chain_init) {
                tree_init->GetEntry(ite_entries);
                timestamp_final = timestamp_init;
                pps_cpt_final = pps_cpt_init;
                pps_info_final = pps_info_init;
                event_cpt_final = event_cpt_init;
                coincidence_final = coincidence_init;
                chip_data_final = chip_data_init;
                analog_trigger_final = analog_trigger_init;
                seu_final = seu_init;
                ch_status_final = ch_status_init;
                ref_channel_final = ref_channel_init;
                if (chain_init == 0) {
                    for (uint64_t i = 0; i<32; i++){
                        sample_p[i] = sample_init[i];
                    }
                    cm_data_final = cm_data_init;
                    tree_init->GetEntry(ite_entries + 1);
                    for (uint64_t i = 0; i<32; i++){
                        sample_n[i] = sample_init[i];
                    }
                } else {
                    for (uint64_t i = 0; i<32; i++){
                        sample_n[i] = sample_init[i];
                    }
                    tree_init->GetEntry(ite_entries + 1);
                    for (uint64_t i = 0; i<32; i++){
                        sample_p[i] = sample_init[i];
                    }
                    cm_data_final = cm_data_init;
                }
                ite_entries += 2;
            } else {
                std::cerr << " == Error == The dataset has 2 values with the same ts, pps and chain. Direct check is needed on iteration " << ite_entries << " and " << ite_entries + 1 << std::endl;
//                tree_init->Show(ite_entries);
//                tree_init->Show(ite_entries + 1);
//                std::exit(EXIT_FAILURE);
                ite_entries += 1;
            }
        } else {
            tree_init->GetEntry(ite_entries);
            timestamp_final = timestamp_init;
            pps_cpt_final = pps_cpt_init;
            pps_info_final = pps_info_init;
            event_cpt_final = event_cpt_init;
            coincidence_final = coincidence_init;
            chip_data_final = chip_data_init;
            analog_trigger_final = analog_trigger_init;
            seu_final = seu_init;
            ch_status_final = ch_status_init;
            ref_channel_final = ref_channel_init;
            cm_data_final = cm_data_init;
            if (chain_init == 0) {
                for (uint64_t i = 0; i<32; i++){
                    sample_p[i] = sample_init[i];
                }
                for (uint64_t i = 0; i<32; i++){
                    sample_n[i] = 0;
                }
            } else {
                for (uint64_t i = 0; i<32; i++){
                    sample_p[i] = 0;
                }
                for (uint64_t i = 0; i<32; i++){
                    sample_n[i] = sample_init[i];
                }
            }
            ite_entries += 1;
        }
        tree_final->Fill();
    }
    if (ite_entries == nentries - 1) {
        tree_init->GetEntry(ite_entries);
        timestamp_final = timestamp_init;
        pps_cpt_final = pps_cpt_init;
        pps_info_final = pps_info_init;
        event_cpt_final = event_cpt_init;
        coincidence_final = coincidence_init;
        chip_data_final = chip_data_init;
        analog_trigger_final = analog_trigger_init;
        seu_final = seu_init;
        ch_status_final = ch_status_init;
        ref_channel_final = ref_channel_init;
        cm_data_final = cm_data_init;
        if (chain_init == 0) {
            for (uint64_t i = 0; i<32; i++){
                sample_p[i] = sample_init[i];
            }
            for (uint64_t i = 0; i<32; i++){
                sample_n[i] = 0;
            }
        } else {
            for (uint64_t i = 0; i<32; i++){
                sample_p[i] = 0;
            }
            for (uint64_t i = 0; i<32; i++){
                sample_n[i] = sample_init[i];
            }
        }
        tree_final->Fill();
    } else {
        std::cout << " == End of DSSSD combination == ";
    }
    file_final->Write();
    file_final->Close();
}

void CorrectTimesIJCLAB(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &glitch_corr_count, uint32_t &minor_corr_count) {
    //
    // Function used to correct the timestamps of the root files.
    //
    // Extraction of the root tree
    TChain *tree_init = new TChain("Events");
    tree_init->Add(init_tree_name.c_str());
//    tree_init->Show(100);

    // Variables used to obtain tree values
    uint32_t ts;
    uint32_t pps_cpt;
    uint32_t gps;
    uint16_t ener_1[32];
    uint16_t ener_2[32];

    // Declaring file_final
    TFile* file_final = nullptr;
    // Declaring tree_final
    TTree* tree_final = nullptr;

    std::string namecalibfile = "Kiruna_data/calib_maud/HG_Piedestaux_&_Corrections_module_VOL_05juin24.txt";
    Double_t maudcorrections[64];
    uint16_t maudpedestals[64];

    // Connect to the branches and opening new root tree according to the detector
    if (detector == "maud") {
        // connect timestamp Maud branches
        tree_init->SetBranchAddress("timestamp",   &ts);
        tree_init->SetBranchAddress("pps_cpt",     &pps_cpt);
        tree_init->SetBranchAddress("pps_info",    &gps);
        // connect energy Maud branches
        tree_init->SetBranchAddress("hg0",   &ener_1);
        tree_init->SetBranchAddress("hg1",   &ener_2);

        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Corrected root file for Maud");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Corrected Maud tree");

        std::ifstream maudcalibfile(namecalibfile); // Open file
        if (!maudcalibfile.is_open()) {
            std::cerr << "Error : Opening of the file failed !" << std::endl;
            return; //
        }

        std::string maudline;
        uint16_t tokencounting;
        size_t start, end;
        size_t channel;
        const std::string delim = " \t"; // " " and tab as delimitor
        std::string corrstring;
        if (std::getline(maudcalibfile, maudline)) {
            std::cout << "Reading Maud calibration file" << std::endl;
            std::cout << "First line giving tokens : " << maudline << std::endl;
        } else {
            std::cerr << "Unable to get the first line of the maud calibration file" << std::endl;
        }
        while (std::getline(maudcalibfile, maudline)) {
            tokencounting = 0;
            start = 0;
            end = 0;

            while ((end = maudline.find_first_of(delim, start)) != std::string::npos) {
                if (end != start) { // Ignore consecutive delimitors
                    if (tokencounting == 0) {
                        channel = std::stoul(maudline.substr(start, end - start));
                    } else if (tokencounting == 4) {
                        corrstring = maudline.substr(start, end - start);
                        std::replace(corrstring.begin(), corrstring.end(), ',', '.');
                        maudcorrections[channel] = std::stod(corrstring);
                    }
                    tokencounting += 1;
                }
                start = end + 1;
            }

            if (start < maudline.size()) {
                maudpedestals[channel] = static_cast<uint16_t>(std::stoi(maudline.substr(start, end - start)));
            }
        }
        maudcalibfile.close();

    } else if (detector == "dssd") {
        // connect timestamp dssd branches
        tree_init->SetBranchAddress("timestamp",   &ts);
        tree_init->SetBranchAddress("pps_cpt",     &pps_cpt);
        tree_init->SetBranchAddress("pps_info",    &gps);
        // connect energy dssd branches
        tree_init->SetBranchAddress("sample_p",   &ener_1);
        tree_init->SetBranchAddress("sample_n",   &ener_2);

        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Corrected root file for DSSD");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Corrected DSSD tree");
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
    //// To validate for energy
    Double_t energy_adc;
    uint16_t adc_channels[64];
    Double_t position[3];
    tree_final->Branch("ts_init", &ts_init, "ts_init/i");
    tree_final->Branch("pps_cpt_init", &pps_cpt_init, "pps_cpt_init/i");
    tree_final->Branch("ts_corr", &ts_corr, "ts_corr/i");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr, "pps_cpt_corr/i");
    tree_final->Branch("gps_corr", &gps_corr, "gps_corr/i");
    tree_final->Branch("time_corr", &time_corr, "time_corr/D");
    //// To validate for energy
    tree_final->Branch("energy_adc", &energy_adc, "energy_adc/D");
    tree_final->Branch("adc_channels", &adc_channels, "adc_channels[64]/s");
    tree_final->Branch("position", &position, "position[3]/D");

    // Variables not saved, used for correcting the time series and creating the new tree
    uint32_t correction_counter = 0;
    uint32_t old_ts=0;
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
        ADCFormatingIJCLAB(ener_1, ener_2, energy_adc, adc_channels, detector, maudcorrections, maudpedestals);
        FindPosition(adc_channels, position, detector);
        tree_final->Fill();
        old_ts = ts;
    }
    file_final->Write();
    file_final->Close();
}

void CorrectTimesCEA(const std::string &init_tree_name, const std::string &final_tree_name, uint32_t &glitch_corr_count, uint32_t &minor_corr_count) {
    //
    // Function used to correct the timestamps of the root files.
    //
    // Extraction of the root tree
    TChain *tree_init = new TChain("D1a");
    tree_init->Add(init_tree_name.c_str());
//    tree_init->Show(100);

    // Variables used to obtain tree values
    Double_t ts;
    Double_t pps_cpt;
    Double_t gps;
    Double_t ener_1[32];
    Double_t ener_2[32];

    // Declaring file_final
    TFile* file_final = nullptr;
    // Declaring tree_final
    TTree* tree_final = nullptr;

    // connect timestamp cea branches
    tree_init->SetBranchAddress("timestamp1",   &ts);
    tree_init->SetBranchAddress("timestamp2",   &pps_cpt);
//        tree_init->SetBranchAddress("pps_info",    &gps);
    gps = 0;
    // connect energy cea branches
    for (uint64_t i = 0; i<32; i++){
        tree_init->SetBranchAddress((std::string("pix") + std::to_string(i + 1)).c_str(),   &ener_1[i]);
        tree_init->SetBranchAddress((std::string("pix") + std::to_string(i + 33)).c_str(),   &ener_2[i]);
    }

    file_exists(final_tree_name);
    file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Corrected root file for CEA");
    if (!file_final || file_final->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    tree_final = new TTree("Events", "Corrected CEA tree");

    // Variables to create the new tree + linking them to the tree
    uint32_t ts_init;
    uint32_t pps_cpt_init;
    uint32_t ts_corr;
    uint32_t pps_cpt_corr;
    uint32_t gps_corr;
    Double_t time_corr;
    //// To validate for energy
    Double_t energy_adc;
    uint16_t adc_channels[64];
    Double_t position[3];
    tree_final->Branch("ts_init", &ts_init, "ts_init/i");
    tree_final->Branch("pps_cpt_init", &pps_cpt_init, "pps_cpt_init/i");
    tree_final->Branch("ts_corr", &ts_corr, "ts_corr/i");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr, "pps_cpt_corr/i");
    tree_final->Branch("gps_corr", &gps_corr, "gps_corr/i");
    tree_final->Branch("time_corr", &time_corr, "time_corr/D");
    //// To validate for energy
    tree_final->Branch("energy_adc", &energy_adc, "energy_adc/D");
    tree_final->Branch("adc_channels", &adc_channels, "adc_channels[64]/s");
    tree_final->Branch("position", &position, "position[3]/D");

    // Variables not saved, used for correcting the time series and creating the new tree
    uint32_t correction_counter = 0;
    uint32_t old_ts=0;
    uint32_t next_ts;
    std::cout << "\n=======================================================================================" << std::endl;
    std::cout << " Correction of cea timestamps "<< std::endl;
    std::cout << "=======================================================================================" << std::endl;
    uint64_t nentries = tree_init->GetEntries();
    for (uint64_t i = 0; i<nentries; i++) {
        // Loading the tree entry
        tree_init->GetEntry(i);
        // Getting the initial, uncorrected, tree values
        pps_cpt_init = static_cast<uint32_t>(pps_cpt);
        ts_init = static_cast<uint32_t>(ts);
        // Starting the correction if needed
        if (ts >= 40.e6) {
            ts_corr = static_cast<uint32_t>(ts) - 40.e6;
            tree_init->GetEntry(i+1);
            next_ts = ts;
            tree_init->GetEntry(i);
            if (old_ts < 40.e6 && next_ts >= 40.e6) {
                correction_counter += 1;
                pps_cpt_corr = static_cast<uint32_t>(pps_cpt) + correction_counter;
                std::cout << "cea correction on entry " << i << ", pps count " << pps_cpt << " to " << pps_cpt_corr << std::endl;
                glitch_corr_count += 1;
            } else if (old_ts < 40.e6 && next_ts < 40.e6){
                std::cout << "ts threshold exceeded but might be due to delay in saving the value, correction applied for this value only, verification is needed to see if there is no shift compared to the GPS clock" << std::endl;
                pps_cpt_corr = static_cast<uint32_t>(pps_cpt) + correction_counter + 1;
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
                pps_cpt_corr = static_cast<uint32_t>(pps_cpt) + correction_counter;
            }
        } else if (ts >= 80.e6) {
            throw std::runtime_error("ERROR : 2 pps incrementations not done in a row, code update needed to solve this");
        } else {
            ts_corr = static_cast<uint32_t>(ts);
            pps_cpt_corr = static_cast<uint32_t>(pps_cpt) + correction_counter;
        }
        gps_corr = static_cast<uint32_t>(gps);
        time_corr = pps_cpt_corr + ts_corr/40.e6;
        ADCFormatingCEA(ener_1, ener_2, energy_adc, adc_channels);
        FindPosition(adc_channels, position, "cea");

        tree_final->Fill();
        old_ts = ts;
    }
    file_final->Write();
    file_final->Close();
}

void ChangeTimeOriginIJCLAB(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, int32_t &abs_gps_time_ref) {
    //
    // Function used to correct the align the timestamp with a common time origin
    //

    // Extraction of the root tree
    file_exists(init_tree_name);
    auto file_init = new TFile(init_tree_name.c_str());
    if (!file_init || file_init->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_init = (TTree*) file_init->Get("Events");
//    std::cout << "KEYS : " << file_init->ls() << std::endl;
//    std::cout << "NOMBRE DENTREES : " << tree_init->GetEntries() << std::endl;

    // Variables to read the initial tree
    uint32_t ts_init;
    uint32_t pps_cpt_init;
    uint32_t ts_corr;
    uint32_t pps_cpt_corr;
    uint32_t gps_corr;
    Double_t time_corr;
    Double_t energy_adc_init;
    uint16_t adc_channels_init[64];
    Double_t position_init[3];

    tree_init->SetBranchAddress("ts_init",       &ts_init);
    tree_init->SetBranchAddress("pps_cpt_init",  &pps_cpt_init);
    tree_init->SetBranchAddress("ts_corr",       &ts_corr);
    tree_init->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr);
    tree_init->SetBranchAddress("gps_corr",      &gps_corr);
    tree_init->SetBranchAddress("time_corr",     &time_corr);
    tree_init->SetBranchAddress("energy_adc",    &energy_adc_init);
    tree_init->SetBranchAddress("adc_channels",  &adc_channels_init);
    tree_init->SetBranchAddress("position",      &position_init);

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    // Opening new root tree according to the detector
    if (detector == "maud") {
//        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for Maud");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time Maud tree");
    } else if (detector == "dssd") {
//        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for DSSD");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time DSSD tree");
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud or dssd" << std::endl;
    }

    // Variables to create the new tree + linking them to the tree
    uint32_t ts_init_noabs;
    int32_t pps_cpt_init_noabs;
    int32_t pps_cpt_corr_noabs;
    int32_t gps_corr_noabs;
    uint32_t ts_corr_abs;
    int32_t pps_cpt_corr_abs;
    int32_t pps_cpt_nocorr_abs;
    int32_t gps_corr_abs;
    Double_t time_corr_abs;
    Double_t energy_adc_final;
    uint16_t adc_channels_final[64];
    Double_t position_final[3];
    Double_t energy_final;
    Double_t temp_final;
    int8_t det_id;

    tree_final->Branch("ts_init", &ts_init_noabs, "ts_init/i");
    tree_final->Branch("pps_cpt_init", &pps_cpt_init_noabs, "pps_cpt_init/I");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr_noabs, "pps_cpt_corr/I");
    tree_final->Branch("gps_corr", &gps_corr_noabs, "gps_corr/I");
    tree_final->Branch("ts_corr_abs", &ts_corr_abs, "ts_corr_abs/i");
    tree_final->Branch("pps_cpt_corr_abs", &pps_cpt_corr_abs, "pps_cpt_corr_abs/I");
    tree_final->Branch("pps_cpt_init_abs", &pps_cpt_nocorr_abs, "pps_cpt_init_abs/I");
    tree_final->Branch("gps_corr_abs", &gps_corr_abs, "gps_corr_abs/I");
    tree_final->Branch("time_corr_abs", &time_corr_abs, "time_corr_abs/D");
    tree_final->Branch("energy_adc", &energy_adc_final, "energy_adc/D");
    tree_final->Branch("adc_channels", &adc_channels_final, "adc_channels[64]/s");
    tree_final->Branch("position", &position_final, "position[3]/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("temperature", &temp_final, "temperature/D");
    tree_final->Branch("det_id", &det_id, "det_id/b");


//    Finds the first entry where the time gps is 1719176452 (97252 in absolute timestamp) (number of second at 23rd june 23:00:52 with origin being 22nd june 20:00:00) return this entry
//    This correction can also be done using the pps_cpt_corr at the first value of 78721 and for pps_cpt_init at the value 78699
//    --  These values of pps were obtained using the GPS time, to adapt this function to the cea case where there is no GPS timestamp
    int32_t ref1 = 1719176452;
    int32_t ref2 = 1719068039;
    uint32_t inc_index = tree_init->GetEntries();
    uint32_t inc_index2 = 0;
    for (uint32_t i=0; i<tree_init->GetEntries(); i++) {
        tree_init->GetEntry(i);
        if (gps_corr == ref1) {
            std::cout << "GPS absolute time of " << ref1 << " incremented at entry " << i << " - Returning this value." << std::endl;
            std::cout << "Absolute time chosen for aligning PPS and GPS : " << gps_corr << std::endl;
            inc_index = i;
            break;
        }
    }
    if (inc_index == tree_init->GetEntries()){
        std::cerr << "ERROR : No entry found for the searched GPS absolute time." << std::endl;
    }

    tree_init->GetEntry(inc_index);
    int32_t abs_pps_time_ref = static_cast<int32_t>(pps_cpt_corr) - gps_corr + abs_gps_time_ref;
    int32_t abs_pps_time_ref2 = 0;
    std::cout << "      Time correction applied to PSS  :  " << abs_pps_time_ref << "      Corresponding entry index : " << inc_index << std::endl;
    std::cout << "      Values of the PPS and GPS value  " << static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref << "    " << gps_corr - abs_gps_time_ref << std::endl;

    if (detector == "maud") {
        for (uint32_t i=0; i<tree_init->GetEntries(); i++) {
            tree_init->GetEntry(i);
            if (gps_corr == ref2) {
                std::cout << "Additional correction value for the first part of the maud data during the ascent phase" << std::endl;
                std::cout << "GPS absolute time of " << ref2 << " incremented at entry " << i << " - Returning this value." << std::endl;
                std::cout << "Absolute time chosen for aligning PPS and GPS : " << gps_corr << std::endl;
                inc_index2 = i;
                break;
            }
        }
        if (inc_index2 == tree_init->GetEntries()){
            std::cerr << "ERROR : No entry found for the SECOND searched GPS absolute time." << std::endl;
        }
        tree_init->GetEntry(inc_index2);
        abs_pps_time_ref2 = static_cast<int32_t>(pps_cpt_corr) - gps_corr + abs_gps_time_ref;
        std::cout << "      Time correction applied to first PSS  :  " << abs_pps_time_ref2 << "      Corresponding entry index : " << inc_index2 << std::endl;
        std::cout << "      Values of the PPS and GPS value  " << static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref2 << "    " << gps_corr - abs_gps_time_ref << std::endl;
    }

    std::vector<int32_t> tstemp_vec;
    std::vector<Double_t> temperature_vec;
    uint64_t temp_vec_init_idx = 0;

    std::vector<int32_t> tscorr_vec;
    std::vector<Double_t> correction_vec;
    uint64_t corr_vec_init_idx = 0;
    Double_t temp_corr;

    Double_t coeffs[64][8];

    if (detector == "maud") {
        std::string tstemp_filename = "Kiruna_data/calib_maud/maud_time_vs_temp_flight.txt";
        ExtractTemperatureValues(temperature_vec, tstemp_vec, tstemp_filename);

        std::string tscorr_filename = "Kiruna_data/calib_maud/fine_calib.txt";
        ExtractCorrectionValues(correction_vec, tscorr_vec, tscorr_filename);

        uint64_t nentries = tree_init->GetEntries();
        for (uint64_t i = 0; i<nentries; i++) {
            tree_init->GetEntry(i);

            ts_init_noabs = ts_init;
            pps_cpt_init_noabs = static_cast<int32_t>(pps_cpt_init);
            pps_cpt_corr_noabs = static_cast<int32_t>(pps_cpt_corr);
            gps_corr_noabs = static_cast<int32_t>(gps_corr);
            ts_corr_abs = ts_corr;
            if (i < 3318949) {
                pps_cpt_corr_abs = static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref2;
                pps_cpt_nocorr_abs = static_cast<int32_t>(pps_cpt_init) - abs_pps_time_ref2;
            } else {
                pps_cpt_corr_abs = static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref;
                pps_cpt_nocorr_abs = static_cast<int32_t>(pps_cpt_init) - abs_pps_time_ref;
            }
            gps_corr_abs = static_cast<int32_t>(gps_corr) - abs_gps_time_ref;
            time_corr_abs = pps_cpt_corr_abs + ts_corr_abs/40.e6;
            energy_adc_final = energy_adc_init;
            for (uint64_t i = 0; i<64; i++){
                adc_channels_final[i] = adc_channels_init[i];
            }
            position_final[0] = position_init[0];
            position_final[1] = position_init[1];
            position_final[2] = position_init[2];

            // Calibration with temperature only
//            FindTemperature(pps_cpt_corr_abs, temp_final, temp_vec_init_idx, tstemp_vec, temperature_vec);
//            temp_corr = 1.;
            // Temperature is set to 0 as there is no use for it in the calibration, so it is not considered to shift the calibration
//            temp_final = 0.;
//            temp_corr = 1.;
            // Calibration with temperature below 19250 and temp corr by fitting over 19250
            // Correction to compensate the effect of temperature that is not uniform
            // A first run is made set to one, to estimate the value of the correction, then the correction is set using the correction file.
            if (pps_cpt_corr_abs < 19250) {
                FindTemperature(pps_cpt_corr_abs, temp_final, temp_vec_init_idx, tstemp_vec, temperature_vec);
                temp_corr = 1.;
            } else {
                temp_final = 0.;
                FindCorrection(pps_cpt_corr_abs, temp_corr, corr_vec_init_idx, tscorr_vec, correction_vec);
            }

            energy_final = ADCtoEnergyMaud(energy_adc_init, temp_final, temp_corr);
            det_id = 7;

            tree_final->Fill();
        }
    } else if (detector == "dssd") {
        std::string tstemp_filename = "Kiruna_data/calib_dssd/dssd_time_vs_temp_flight.txt";
        ExtractTemperatureValues(temperature_vec, tstemp_vec, tstemp_filename);

        std::string dssd_p_coefs_filename = "Kiruna_data/calib_dssd/pside_v3.calib";
        std::string dssd_n_coefs_filename = "Kiruna_data/calib_dssd/nside_v3.calib";
        ExtractDSSDCoefs(coeffs, dssd_p_coefs_filename, dssd_n_coefs_filename);

        Double_t adc_channels_calib[64];
        Double_t adc_channels_calib_corr[64];
//        Double_t adc_channels_init_double[64];
        tree_final->Branch("adc_channels_calib", &adc_channels_calib, "adc_channels_calib[64]/D");
//        tree_final->Branch("adc_channels_cali_corr", &adc_channels_corr, "adc_channels_cali_corr[64]/D");
        tree_final->Branch("adc_channels_calib_corr", &adc_channels_calib_corr, "adc_channels_calib_corr[64]/D");

        uint64_t nentries = tree_init->GetEntries();
        for (uint64_t i = 0; i<nentries; i++) {
            tree_init->GetEntry(i);

            ts_init_noabs = ts_init;
            pps_cpt_init_noabs = static_cast<int32_t>(pps_cpt_init);
            pps_cpt_corr_noabs = static_cast<int32_t>(pps_cpt_corr);
            gps_corr_noabs = static_cast<int32_t>(gps_corr);
            ts_corr_abs = ts_corr;
            pps_cpt_corr_abs = static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref;
            pps_cpt_nocorr_abs = static_cast<int32_t>(pps_cpt_init) - abs_pps_time_ref;
            gps_corr_abs = static_cast<int32_t>(gps_corr) - abs_gps_time_ref;
            time_corr_abs = pps_cpt_corr_abs + ts_corr_abs/40.e6;
            energy_adc_final = energy_adc_init;
            for (uint64_t i = 0; i<64; i++){
                adc_channels_final[i] = adc_channels_init[i];
            }
            position_final[0] = position_init[0];
            position_final[1] = position_init[1];
            position_final[2] = position_init[2];

            FindTemperature(pps_cpt_corr_abs, temp_final, temp_vec_init_idx, tstemp_vec, temperature_vec);
//            temp_final = 0.;
//            ADCTempCorrectionDSSD(adc_channels_init, adc_channels_corr, temp_final);
            // This line for applying the temperature shift before calibration
            energy_final = ADCtoEnergyDSSD(adc_channels_init, adc_channels_calib, adc_channels_calib_corr, temp_final, coeffs); //!! c'est comme si c'était à 25°
            // This line for applying the calibration to raw channels
//            for (uint64_t i = 0; i<64; i++){
//                adc_channels_init_double[i] = static_cast<Double_t>(adc_channels_init[i]);
//            }
//            energy_final = ADCtoEnergyDSSD(adc_channels_init_double, adc_channels_calib, temp_final, coeffs);
            det_id = 2;

            tree_final->Fill();
        }
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud or dssd" << std::endl;
    }
    file_init->Close();
    file_final->Write();
    file_final->Close();
}

void ChangeTimeOriginCEA(const std::string &init_tree_name, const std::string &final_tree_name, int32_t &abs_pps_time_ref) {
    //
    // Function used to correct the align the timestamp with a common time origin
    //

    // Extraction of the root tree
    file_exists(init_tree_name);
    auto file_init = new TFile(init_tree_name.c_str());
    if (!file_init || file_init->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_init = (TTree*) file_init->Get("Events");

    // Variables to read the initial tree
    uint32_t ts_init;
    uint32_t pps_cpt_init;
    uint32_t ts_corr;
    uint32_t pps_cpt_corr;
    uint32_t gps_corr;
    Double_t time_corr;
    Double_t energy_adc_init;
    uint16_t adc_channels_init[64];
    Double_t position_init[3];

    tree_init->SetBranchAddress("ts_init",       &ts_init);
    tree_init->SetBranchAddress("pps_cpt_init",  &pps_cpt_init);
    tree_init->SetBranchAddress("ts_corr",       &ts_corr);
    tree_init->SetBranchAddress("pps_cpt_corr",  &pps_cpt_corr);
    tree_init->SetBranchAddress("gps_corr",      &gps_corr);
    tree_init->SetBranchAddress("time_corr",     &time_corr);
    tree_init->SetBranchAddress("energy_adc",    &energy_adc_init);
    tree_init->SetBranchAddress("adc_channels",  &adc_channels_init);
    tree_init->SetBranchAddress("position",      &position_init);

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    // Opening new root tree
    file_exists(final_tree_name);
    file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for CEA");
    if (!file_final || file_final->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    tree_final = new TTree("Events", "Absolute time CEA tree");

    // Variables to create the new tree + linking them to the tree
    uint32_t ts_init_noabs;
    int32_t pps_cpt_init_noabs;
    int32_t pps_cpt_corr_noabs;
    int32_t gps_corr_noabs;
    uint32_t ts_corr_abs;
    int32_t pps_cpt_corr_abs;
    int32_t pps_cpt_nocorr_abs;
    int32_t gps_corr_abs;
    Double_t time_corr_abs;
    Double_t energy_adc_final;
    uint16_t adc_channels_final[64];
    Double_t position_final[3];
    Double_t energy_final;
    Double_t temp_final;
    int8_t det_id;

    Double_t adc_channels_calib[64];
    Double_t adc_channels_corr[64];
//    Double_t adc_channels_init_double[64];

    tree_final->Branch("ts_init", &ts_init_noabs, "ts_init/i");
    tree_final->Branch("pps_cpt_init", &pps_cpt_init_noabs, "pps_cpt_init/I");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr_noabs, "pps_cpt_corr/I");
    tree_final->Branch("gps_corr", &gps_corr_noabs, "gps_corr/I");
    tree_final->Branch("ts_corr_abs", &ts_corr_abs, "ts_corr_abs/i");
    tree_final->Branch("pps_cpt_corr_abs", &pps_cpt_corr_abs, "pps_cpt_corr_abs/I");
    tree_final->Branch("pps_cpt_init_abs", &pps_cpt_nocorr_abs, "pps_cpt_init_abs/I");
    tree_final->Branch("gps_corr_abs", &gps_corr_abs, "gps_corr_abs/I");
    tree_final->Branch("time_corr_abs", &time_corr_abs, "time_corr_abs/D");
    tree_final->Branch("energy_adc", &energy_adc_final, "energy_adc/D");
    tree_final->Branch("adc_channels", &adc_channels_final, "adc_channels[64]/s");
    tree_final->Branch("position", &position_final, "position[3]/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("temperature", &temp_final, "temperature/D");
    tree_final->Branch("det_id", &det_id, "det_id/b");

    tree_final->Branch("adc_channels_calib", &adc_channels_calib, "adc_channels_calib[64]/D");
    tree_final->Branch("adc_channels_corr", &adc_channels_corr, "adc_channels_corr[64]/D");

//    Finds the first entry where the absolute time gps is 97252 (number of second at 23rd june 23:00:52 with origin being 22nd june 20:00:00) return this entry
//    This correction can also be done using the pps_cpt_corr at the first value of 78721 and for pps_cpt_init at the value 78699
//    --  These values of pps were obtained using the GPS time, to adapt this function to the cea case where there is no GPS timestamp
//    uint32_t ref = 97252;
//    uint32_t inc_index = tree_init->GetEntries();
//    for (uint32_t i=0; i<tree_init->GetEntries(); i++) {
//        tree_init->GetEntry(i);
//        if (gps_corr - abs_gps_time_ref == ref){
//            std::cout << "GPS absolute time of " << ref << " incremented at entry " << i << " - Returning this value." << std::endl;
//            std::cout << "Absolute time chosen for aligning PPS and GPS : " << gps_corr - abs_gps_time_ref << std::endl;
//            inc_index = i;
//            break;
//            }
//    }
//    if (inc_index == tree_init->GetEntries()){
//        std::cerr << "ERROR : No entry found for the searched GPS absolute time." << std::endl;
//    }

//    tree_init->GetEntry(inc_index);
//    uint32_t abs_pps_time_ref = pps_cpt_corr - gps_corr + abs_gps_time_ref;
//    std::cout << "      Time correction applied to PSS  :  " << abs_pps_time_ref << "      Corresponding entry index : " << inc_index << std::endl;
//    std::cout << "      Values of the PPS and GPS value  " << pps_cpt_corr - abs_pps_time_ref << "    " << gps_corr - abs_gps_time_ref << std::endl;
    std::cout << "      Time correction applied to PSS  :  " << abs_pps_time_ref << std::endl;

//    std::vector<uint32_t> tstemp_vec;
//    std::vector<Double_t> temperature_vec;
//    uint64_t temp_vec_init_idx = 0;
//    std::string tstemp_filename = "Kiruna_data/calib_cea/cea_time_vs_temp_flight.txt";
//    ExtractTemperatureValues(temperature_vec, tstemp_vec, tstemp_filename);

    Double_t coeffs[64][8];

    std::string cea_p_coefs_filename = "Kiruna_data/calib_cea/pside.calib";
    std::string cea_n_coefs_filename = "Kiruna_data/calib_cea/nside.calib";
    ExtractDSSDCoefs(coeffs, cea_p_coefs_filename, cea_n_coefs_filename);

    uint64_t nentries = tree_init->GetEntries();
    for (uint64_t i = 0; i<nentries; i++) {
        tree_init->GetEntry(i);

        ts_init_noabs = ts_init;
        pps_cpt_init_noabs = static_cast<int32_t>(pps_cpt_init);
        pps_cpt_corr_noabs = static_cast<int32_t>(pps_cpt_corr);
        gps_corr_noabs = static_cast<int32_t>(gps_corr);
        ts_corr_abs = ts_corr;
        pps_cpt_corr_abs = static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref;
        pps_cpt_nocorr_abs = static_cast<int32_t>(pps_cpt_init) - abs_pps_time_ref;
        gps_corr_abs = pps_cpt_corr_abs;
        time_corr_abs = pps_cpt_corr_abs + ts_corr_abs/40.e6;
        energy_adc_final = energy_adc_init;
        for (uint64_t i = 0; i<64; i++){
            adc_channels_final[i] = adc_channels_init[i];
        }
        position_final[0] = position_init[0];
        position_final[1] = position_init[1];
        position_final[2] = position_init[2];

//        FindTemperature(pps_cpt_corr_abs, temp_final, temp_vec_init_idx, tstemp_vec, temperature_vec);
        temp_final = 0.;
        ADCTempCorrectionCEA(adc_channels_init, adc_channels_corr, temp_final);
        // This line for applying the temperature shift before calibration
        energy_final = ADCtoEnergyCEA(adc_channels_init, adc_channels_calib, 25, coeffs); //!! c'est comme si c'était à 25°
        // This line for applying the calibration to raw channels
//            for (uint64_t i = 0; i<64; i++){
//                adc_channels_init_double[i] = static_cast<Double_t>(adc_channels_init[i]);
//            }
//            energy_final = ADCtoEnergyDSSD(adc_channels_init_double, adc_channels_calib, temp_final, coeffs);
        det_id = 1;
        tree_final->Fill();
    }
    file_init->Close();
    file_final->Write();
    file_final->Close();
}

void ChangeTimeOriginUCD(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, int32_t &abs_gps_time_ref) {
    //
    // Function used to correct the align the timestamp with a common time origin
    //

    // Extraction of the root tree
    file_exists(init_tree_name);
    auto file_init = new TFile(init_tree_name.c_str());
    if (!file_init || file_init->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_init = (TTree*) file_init->Get("Events");
//    std::cout << "KEYS : " << file_init->ls() << std::endl;
//    std::cout << "NOMBRE DENTREES : " << tree_init->GetEntries() << std::endl;

    // Variables to read the initial tree
    int32_t pps_cpt_init;
    Long64_t pps_cpt_corr;
    int32_t gps_corr;
    Double_t time_corr;
    Double_t energy_init;
    Short_t position_init;
    tree_init->SetBranchAddress("Time_sec",  &pps_cpt_init);
    tree_init->SetBranchAddress("Time_sec_corrected",  &pps_cpt_corr);
    tree_init->SetBranchAddress("Time_gps",      &gps_corr);
    tree_init->SetBranchAddress("Time",     &time_corr);
    tree_init->SetBranchAddress("Energy",        &energy_init);
    tree_init->SetBranchAddress("Argmax",        &position_init);

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    // Opening new root tree according to the detector
    if (detector == "ucda") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for UCDA");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time UCDA tree");
    } else if (detector == "ucdb") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for UCDB");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time UCDB tree");
    } else if (detector == "ucdc") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for UCDC");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time UCDC tree");
    } else if (detector == "ucdd") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for UCDD");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Absolute time UCDD tree");
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be ucda, ucdb, ucdc or ucdd" << std::endl;
    }

    // Variables to create the new tree + linking them to the tree
    int32_t pps_cpt_init_noabs;
    int32_t pps_cpt_corr_noabs;
    int32_t gps_corr_noabs;
    int32_t pps_cpt_corr_abs;
    int32_t pps_cpt_nocorr_abs;
    int32_t gps_corr_abs;
    Double_t time_corr_abs;
    Double_t energy_final;
    Double_t position_final[3];
    tree_final->Branch("pps_cpt_init", &pps_cpt_init_noabs, "pps_cpt_init/I");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr_noabs, "pps_cpt_corr/I");
    tree_final->Branch("gps_corr", &gps_corr_noabs, "gps_corr/I");
    tree_final->Branch("pps_cpt_corr_abs", &pps_cpt_corr_abs, "pps_cpt_corr_abs/I");
    tree_final->Branch("pps_cpt_init_abs", &pps_cpt_nocorr_abs, "pps_cpt_init_abs/I");
    tree_final->Branch("gps_corr_abs", &gps_corr_abs, "gps_corr_abs/I");
    tree_final->Branch("time_corr_abs", &time_corr_abs, "time_corr_abs/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("position", &position_final, "position[3]/D");

    //    Finds the first entry where the time gps is 1719176452 (97252 in absolute timestamp) (number of second at 23rd june 23:00:52 with origin being 22nd june 20:00:00) return this entry
    int32_t ref = 1719176452;
    uint32_t inc_index = tree_init->GetEntries();
    for (uint32_t i=0; i<tree_init->GetEntries(); i++) {
        tree_init->GetEntry(i);
        if (gps_corr == ref){
            std::cout << "GPS absolute time of " << ref << " incremented at entry " << i << " - Returning this value." << std::endl;
            std::cout << "Absolute time chosen for aligning PPS and GPS : " << gps_corr << std::endl;
            inc_index = i;
            break;
            }
    }
    if (inc_index == tree_init->GetEntries()){
        std::cerr << "ERROR : No entry found for the searched GPS absolute time." << std::endl;
    }

    tree_init->GetEntry(inc_index);
    int32_t abs_pps_time_ref = static_cast<int32_t>(pps_cpt_corr) - gps_corr + abs_gps_time_ref;
    std::cout << "      Time correction applied to PSS  :  " << abs_pps_time_ref << "      Corresponding entry index : " << inc_index << std::endl;
    std::cout << "      Values of the PPS and GPS value  " << static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref << "    " << gps_corr - abs_gps_time_ref << std::endl;

    uint64_t nentries = tree_init->GetEntries();
    for (uint64_t i = 0; i<nentries; i++) {
        tree_init->GetEntry(i);

        pps_cpt_init_noabs = pps_cpt_init;
        pps_cpt_corr_noabs = static_cast<int32_t>(pps_cpt_corr);
        gps_corr_noabs = gps_corr;
        pps_cpt_corr_abs = static_cast<int32_t>(pps_cpt_corr) - abs_pps_time_ref;
        pps_cpt_nocorr_abs = pps_cpt_init - abs_pps_time_ref;
        gps_corr_abs = gps_corr - abs_gps_time_ref;
        time_corr_abs = time_corr;
        energy_final = energy_init;
        GetUCDPosition(position_final, position_init);
        tree_final->Fill();
    }
    file_init->Close();
    file_final->Write();
    file_final->Close();
}

void BuildFinalTreeIJCLAB(const std::string &init_tree_name, const std::string &final_tree_name, const std::vector<uint64_t> coinc_idx_ijclab, const std::vector<uint64_t> coinc_idx_cea, const std::vector<uint64_t> coinc_idx_ucda, const std::vector<uint64_t> coinc_idx_ucdb, const std::vector<uint64_t> coinc_idx_ucdc, const std::vector<uint64_t> coinc_idx_ucdd, std::string detector) {
    //
    // Function used to add coincidences flags in the tree
    //

    // Extraction of the root tree
    file_exists(init_tree_name);
    auto file_init = new TFile(init_tree_name.c_str());
    if (!file_init || file_init->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_init = (TTree*) file_init->Get("Events");

    // Variables to read the initial tree
    Double_t time_init;
    Double_t position_init[3];
    Double_t energy_init;

    tree_init->SetBranchAddress("time_corr_abs", &time_init);
    tree_init->SetBranchAddress("position", &position_init);
    tree_init->SetBranchAddress("energy", &energy_init);

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    // Opening new root tree according to the detector
    if (detector == "maud") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Final root file for Maud");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Final Maud tree");
    } else if (detector == "dssd") {
        file_exists(final_tree_name);
        file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Final root file for DSSD");
        if (!file_final || file_final->IsZombie()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
        tree_final = new TTree("Events", "Final DSSD tree");
    } else {
        std::cerr << " ERROR : Unknown detector : " << detector << " \n detector must be maud or dssd" << std::endl;
    }

    // Variables to create the new tree + linking them to the tree
    Double_t time_final;
    Double_t position_final[3];
    Double_t energy_final;
    uint8_t coinc_ijclab;
    uint8_t coinc_cea;
    uint8_t coinc_ucda;
    uint8_t coinc_ucdb;
    uint8_t coinc_ucdc;
    uint8_t coinc_ucdd;
    uint8_t chain_len;

    tree_final->Branch("time", &time_final, "time/D");
    tree_final->Branch("position", &position_final, "position[3]/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("coinc_ijclab", &coinc_ijclab, "coinc_ijclab/b");
    tree_final->Branch("coinc_cea", &coinc_cea, "coinc_cea/b");
    tree_final->Branch("coinc_ucda", &coinc_ucda, "coinc_ucda/b");
    tree_final->Branch("coinc_ucdb", &coinc_ucdb, "coinc_ucdb/b");
    tree_final->Branch("coinc_ucdc", &coinc_ucdc, "coinc_ucdc/b");
    tree_final->Branch("coinc_ucdd", &coinc_ucdd, "coinc_ucdd/b");
    tree_final->Branch("chain_len", &chain_len, "chain_len/b");

    uint64_t nentries = tree_init->GetEntries();
    // Using vectors to easily fill the tree
    std::vector<uint8_t> vec_coinc_ijclab(nentries, 0);
    std::vector<uint8_t> vec_coinc_cea(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucda(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucdb(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucdc(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucdd(nentries, 0);

    for (uint64_t i = 0; i<coinc_idx_ijclab.size(); i++) {
        vec_coinc_ijclab[coinc_idx_ijclab[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_cea.size(); i++) {
        vec_coinc_cea[coinc_idx_cea[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucda.size(); i++) {
        vec_coinc_ucda[coinc_idx_ucda[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucdb.size(); i++) {
        vec_coinc_ucdb[coinc_idx_ucdb[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucdc.size(); i++) {
        vec_coinc_ucdc[coinc_idx_ucdc[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucdd.size(); i++) {
        vec_coinc_ucdd[coinc_idx_ucdd[i]] = 1;
    }

    for (uint64_t i = 0; i<nentries; i++) {
        tree_init->GetEntry(i);

        time_final = time_init;
        position_final[0] = position_init[0];
        position_final[1] = position_init[1];
        position_final[2] = position_init[2];
        energy_final = energy_init;

        coinc_ijclab = vec_coinc_ijclab[i];
        coinc_cea = vec_coinc_cea[i];
        coinc_ucda = vec_coinc_ucda[i];
        coinc_ucdb = vec_coinc_ucdb[i];
        coinc_ucdc = vec_coinc_ucdc[i];
        coinc_ucdd = vec_coinc_ucdd[i];
        chain_len = vec_coinc_ijclab[i] + vec_coinc_cea[i] + vec_coinc_ucda[i] + vec_coinc_ucdb[i] + vec_coinc_ucdc[i] + vec_coinc_ucdd[i];
        tree_final->Fill();
    }

    file_init->Close();
    file_final->Write();
    file_final->Close();
}

void BuildFinalTreeUCD(const std::string &tree_name_a, const std::string &tree_name_b, const std::string &tree_name_c, const std::string &tree_name_d, const std::string &final_tree_name, const std::vector<uint64_t> coinc_idx_ijclab1, const std::vector<uint64_t> coinc_idx_ijclab2, const std::vector<uint64_t> coinc_idx_cea, const std::vector<uint64_t> coinc_idx_a_b, const std::vector<uint64_t> coinc_idx_a_c, const std::vector<uint64_t> coinc_idx_a_d, const std::vector<uint64_t> coinc_idx_b_a, const std::vector<uint64_t> coinc_idx_b_c, const std::vector<uint64_t> coinc_idx_b_d, const std::vector<uint64_t> coinc_idx_c_a, const std::vector<uint64_t> coinc_idx_c_b, const std::vector<uint64_t> coinc_idx_c_d, const std::vector<uint64_t> coinc_idx_d_a, const std::vector<uint64_t> coinc_idx_d_b, const std::vector<uint64_t> coinc_idx_d_c) {
    // Loading tree a
    file_exists(tree_name_a);
    auto file_a = new TFile(tree_name_a.c_str());
    if (!file_a || file_a->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_a = (TTree*) file_a->Get("Events");

    // Variables to read the initial tree
    Double_t time_corr_abs_a;
    Double_t position_a[3];
    Double_t energy_final_a;
    tree_a->SetBranchAddress("time_corr_abs", &time_corr_abs_a);
    tree_a->SetBranchAddress("position", &position_a);
    tree_a->SetBranchAddress("energy", &energy_final_a);

    // Loading tree b
    file_exists(tree_name_b);
    auto file_b = new TFile(tree_name_b.c_str());
    if (!file_b || file_b->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_b = (TTree*) file_b->Get("Events");

    // Variables to read the initial tree
    Double_t time_corr_abs_b;
    Double_t position_b[3];
    Double_t energy_final_b;
    tree_b->SetBranchAddress("time_corr_abs", &time_corr_abs_b);
    tree_b->SetBranchAddress("position", &position_b);
    tree_b->SetBranchAddress("energy", &energy_final_b);

    // Loading tree c
    file_exists(tree_name_c);
    auto file_c = new TFile(tree_name_c.c_str());
    if (!file_c || file_c->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_c = (TTree*) file_c->Get("Events");

    // Variables to read the initial tree
    Double_t time_corr_abs_c;
    Double_t position_c[3];
    Double_t energy_final_c;
    tree_c->SetBranchAddress("time_corr_abs", &time_corr_abs_c);
    tree_c->SetBranchAddress("position", &position_c);
    tree_c->SetBranchAddress("energy", &energy_final_c);

    // Loading tree d
    file_exists(tree_name_d);
    auto file_d = new TFile(tree_name_d.c_str());
    if (!file_d || file_d->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_d = (TTree*) file_d->Get("Events");

    // Variables to read the initial tree
    Double_t time_corr_abs_d;
    Double_t position_d[3];
    Double_t energy_final_d;
    tree_d->SetBranchAddress("time_corr_abs", &time_corr_abs_d);
    tree_d->SetBranchAddress("position", &position_d);
    tree_d->SetBranchAddress("energy", &energy_final_d);


    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for UCD");
    if (!file_final || file_final->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    tree_final = new TTree("Events", "Absolute time UCD tree");

    // Variables of the new tree
    Double_t time_final;
    Double_t position_final[3];
    Double_t energy_final;
    uint8_t det_id;
    uint8_t coinc_ijclab1;
    uint8_t coinc_ijclab2;
    uint8_t coinc_cea;
    uint8_t coinc_ucda;
    uint8_t coinc_ucdb;
    uint8_t coinc_ucdc;
    uint8_t coinc_ucdd;
    uint8_t chain_len;

    tree_final->Branch("time", &time_final, "time/D");
    tree_final->Branch("position", &position_final, "position[3]/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("det_id", &det_id, "det_id/b");
    tree_final->Branch("coinc_ijclab1", &coinc_ijclab1, "coinc_ijclab1/b");
    tree_final->Branch("coinc_ijclab2", &coinc_ijclab2, "coinc_ijclab2/b");
    tree_final->Branch("coinc_cea", &coinc_cea, "coinc_cea/b");
    tree_final->Branch("coinc_ucda", &coinc_ucda, "coinc_ucda/b");
    tree_final->Branch("coinc_ucdb", &coinc_ucdb, "coinc_ucdb/b");
    tree_final->Branch("coinc_ucdc", &coinc_ucdc, "coinc_ucdc/b");
    tree_final->Branch("coinc_ucdd", &coinc_ucdd, "coinc_ucdd/b");
    tree_final->Branch("chain_len", &chain_len, "chain_len/b");

    uint64_t nentries_a = tree_a->GetEntries();
    // Using vectors to easily fill the tree a
    std::vector<uint8_t> vec_coinc_ijclab1_a(nentries_a, 0);
    std::vector<uint8_t> vec_coinc_ijclab2_a(nentries_a, 0);
    std::vector<uint8_t> vec_coinc_cea_a(nentries_a, 0);
    std::vector<uint8_t> vec_coinc_a_b(nentries_a, 0);
    std::vector<uint8_t> vec_coinc_a_c(nentries_a, 0);
    std::vector<uint8_t> vec_coinc_a_d(nentries_a, 0);
    for (uint64_t i = 0; i<coinc_idx_ijclab1.size(); i++) {
        vec_coinc_ijclab1_a[coinc_idx_ijclab1[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ijclab2.size(); i++) {
        vec_coinc_ijclab2_a[coinc_idx_ijclab2[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_cea.size(); i++) {
        vec_coinc_cea_a[coinc_idx_cea[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_a_b.size(); i++) {
        vec_coinc_a_b[coinc_idx_a_b[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_a_c.size(); i++) {
        vec_coinc_a_c[coinc_idx_a_c[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_a_d.size(); i++) {
        vec_coinc_a_d[coinc_idx_a_d[i]] = 1;
    }

    uint64_t nentries_b = tree_b->GetEntries();
    // Using vectors to easily fill the tree b
    std::vector<uint8_t> vec_coinc_ijclab1_b(nentries_b, 0);
    std::vector<uint8_t> vec_coinc_ijclab2_b(nentries_b, 0);
    std::vector<uint8_t> vec_coinc_cea_b(nentries_b, 0);
    std::vector<uint8_t> vec_coinc_b_a(nentries_b, 0);
    std::vector<uint8_t> vec_coinc_b_c(nentries_b, 0);
    std::vector<uint8_t> vec_coinc_b_d(nentries_b, 0);
    for (uint64_t i = 0; i<coinc_idx_ijclab1.size(); i++) {
        vec_coinc_ijclab1_b[coinc_idx_ijclab1[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ijclab2.size(); i++) {
        vec_coinc_ijclab2_b[coinc_idx_ijclab2[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_cea.size(); i++) {
        vec_coinc_cea_b[coinc_idx_cea[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_b_a.size(); i++) {
        vec_coinc_b_a[coinc_idx_b_a[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_b_c.size(); i++) {
        vec_coinc_b_c[coinc_idx_b_c[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_b_d.size(); i++) {
        vec_coinc_b_d[coinc_idx_b_d[i]] = 1;
    }

    uint64_t nentries_c = tree_c->GetEntries();
    // Using vectors to easily fill the tree c
    std::vector<uint8_t> vec_coinc_ijclab1_c(nentries_c, 0);
    std::vector<uint8_t> vec_coinc_ijclab2_c(nentries_c, 0);
    std::vector<uint8_t> vec_coinc_cea_c(nentries_c, 0);
    std::vector<uint8_t> vec_coinc_c_a(nentries_c, 0);
    std::vector<uint8_t> vec_coinc_c_b(nentries_c, 0);
    std::vector<uint8_t> vec_coinc_c_d(nentries_c, 0);
    for (uint64_t i = 0; i<coinc_idx_ijclab1.size(); i++) {
        vec_coinc_ijclab1_c[coinc_idx_ijclab1[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ijclab2.size(); i++) {
        vec_coinc_ijclab2_c[coinc_idx_ijclab2[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_cea.size(); i++) {
        vec_coinc_cea_c[coinc_idx_cea[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_c_a.size(); i++) {
        vec_coinc_c_a[coinc_idx_c_a[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_c_b.size(); i++) {
        vec_coinc_c_b[coinc_idx_c_b[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_c_d.size(); i++) {
        vec_coinc_c_d[coinc_idx_c_d[i]] = 1;
    }

    uint64_t nentries_d = tree_d->GetEntries();
    // Using vectors to easily fill the tree d
    std::vector<uint8_t> vec_coinc_ijclab1_d(nentries_d, 0);
    std::vector<uint8_t> vec_coinc_ijclab2_d(nentries_d, 0);
    std::vector<uint8_t> vec_coinc_cea_d(nentries_d, 0);
    std::vector<uint8_t> vec_coinc_d_a(nentries_d, 0);
    std::vector<uint8_t> vec_coinc_d_b(nentries_d, 0);
    std::vector<uint8_t> vec_coinc_d_c(nentries_d, 0);
    for (uint64_t i = 0; i<coinc_idx_ijclab1.size(); i++) {
        vec_coinc_ijclab1_d[coinc_idx_ijclab1[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ijclab2.size(); i++) {
        vec_coinc_ijclab2_d[coinc_idx_ijclab2[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_cea.size(); i++) {
        vec_coinc_cea_d[coinc_idx_cea[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_d_a.size(); i++) {
        vec_coinc_d_a[coinc_idx_d_a[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_d_b.size(); i++) {
        vec_coinc_d_b[coinc_idx_d_b[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_d_c.size(); i++) {
        vec_coinc_d_c[coinc_idx_d_c[i]] = 1;
    }

    // Combining the UCD trees
    std::vector<Double_t> timevec;
    timevec.reserve(nentries_a + nentries_b + nentries_c + nentries_d);

    std::cout << "UCDA number of event : " << nentries_a << std::endl;
    for (uint64_t i = 0; i<nentries_a; i++) {
        tree_a->GetEntry(i);
        timevec.push_back(time_corr_abs_a);
    }

    std::cout << "UCDB number of event : " << nentries_b << std::endl;
    for (uint64_t i = 0; i<nentries_b; i++) {
        tree_b->GetEntry(i);
        timevec.push_back(time_corr_abs_b);
    }

    std::cout << "UCDC number of event : " << nentries_c << std::endl;
    for (uint64_t i = 0; i<nentries_c; i++) {
        tree_c->GetEntry(i);
        timevec.push_back(time_corr_abs_c);
    }

    std::cout << "UCDD number of event : " << nentries_d << std::endl;
    for (uint64_t i = 0; i<nentries_d; i++) {
        tree_d->GetEntry(i);
        timevec.push_back(time_corr_abs_d);
    }

    uint64_t nentries_final = timevec.size();
    std::vector<size_t> sorted_index = sort_index(timevec);
    uint64_t sorted_ite;
    for (uint64_t i = 0; i<sorted_index.size(); i++) {
        sorted_ite = sorted_index[i];
        // Getting the entry in sorted order
        if (sorted_ite < nentries_a) {
            tree_a->GetEntry(sorted_ite);
            time_final = time_corr_abs_a;
            position_final[0] = position_a[0];
            position_final[1] = position_a[1];
            position_final[2] = position_a[2];
            energy_final = energy_final_a;
            det_id = 0;
            coinc_ijclab1 = vec_coinc_ijclab1_a[i];
            coinc_ijclab2 = vec_coinc_ijclab2_a[i];
            coinc_cea = vec_coinc_cea_a[i];
            coinc_ucda = 0;
            coinc_ucdb = vec_coinc_a_b[i];
            coinc_ucdc = vec_coinc_a_c[i];
            coinc_ucdd = vec_coinc_a_d[i];
            chain_len = coinc_ijclab1 + coinc_ijclab2 + coinc_cea + coinc_ucda + coinc_ucdb + coinc_ucdc + coinc_ucdd;
        } else if (sorted_ite < nentries_a + nentries_b) {
            tree_b->GetEntry(sorted_ite - nentries_a);
            time_final = time_corr_abs_b;
            position_final[0] = position_b[0];
            position_final[1] = position_b[1];
            position_final[2] = position_b[2];
            energy_final = energy_final_b;
            det_id = 0;
            coinc_ijclab1 = vec_coinc_ijclab1_b[i];
            coinc_ijclab2 = vec_coinc_ijclab2_b[i];
            coinc_cea = vec_coinc_cea_b[i];
            coinc_ucda = vec_coinc_b_a[i];
            coinc_ucdb = 0;
            coinc_ucdc = vec_coinc_b_c[i];
            coinc_ucdd = vec_coinc_b_d[i];
        } else if (sorted_ite < nentries_a + nentries_b + nentries_c) {
            tree_c->GetEntry(sorted_ite - nentries_a - nentries_b);
            time_final = time_corr_abs_c;
            position_final[0] = position_c[0];
            position_final[1] = position_c[1];
            position_final[2] = position_c[2];
            energy_final = energy_final_c;
            det_id = 0;
            coinc_ijclab1 = vec_coinc_ijclab1_c[i];
            coinc_ijclab2 = vec_coinc_ijclab2_c[i];
            coinc_cea = vec_coinc_cea_c[i];
            coinc_ucda = vec_coinc_c_a[i];
            coinc_ucdb = vec_coinc_c_b[i];
            coinc_ucdc = 0;
            coinc_ucdd = vec_coinc_c_d[i];
        } else {
            tree_d->GetEntry(sorted_ite - nentries_a - nentries_b - nentries_c);
            time_final = time_corr_abs_d;
            position_final[0] = position_d[0];
            position_final[1] = position_d[1];
            position_final[2] = position_d[2];
            energy_final = energy_final_d;
            det_id = 0;
            coinc_ijclab1 = vec_coinc_ijclab1_d[i];
            coinc_ijclab2 = vec_coinc_ijclab2_d[i];
            coinc_cea = vec_coinc_cea_d[i];
            coinc_ucda = vec_coinc_d_a[i];
            coinc_ucdb = vec_coinc_d_b[i];
            coinc_ucdc = vec_coinc_d_c[i];
            coinc_ucdd = 0;
        }

        tree_final->Fill();
        if (i % 100000 == 0) {
            std::cout << "\rOrdering the events in the new tree : entry " << i << " over " << nentries_final << std::flush;
        }
    }

    file_final->Write();

    std::cout << "size time vector vs size final tree " << timevec.size() << " - " << tree_final->GetEntries() << std::endl;

    file_a->Close();
    file_b->Close();
    file_c->Close();
    file_d->Close();
    file_final->Close();

}

void BuildFinalTreeCEA(const std::string &init_tree_name, const std::string &final_tree_name, const std::vector<uint64_t> coinc_idx_ijclab1, const std::vector<uint64_t> coinc_idx_ijclab2, const std::vector<uint64_t> coinc_idx_ucda, const std::vector<uint64_t> coinc_idx_ucdb, const std::vector<uint64_t> coinc_idx_ucdc, const std::vector<uint64_t> coinc_idx_ucdd) {
    //
    // Function used to add coincidences flags in the tree
    //

    // Extraction of the root tree
    file_exists(init_tree_name);
    auto file_init = new TFile(init_tree_name.c_str());
    if (!file_init || file_init->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_init = (TTree*) file_init->Get("Events");

    // Variables to read the initial tree
    Double_t time_init;
    Double_t position_init[3];
    Double_t energy_init;

    tree_init->SetBranchAddress("time_corr_abs", &time_init);
    tree_init->SetBranchAddress("position", &position_init);
    tree_init->SetBranchAddress("energy", &energy_init);

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    // Opening new root tree according to the detector
    file_exists(final_tree_name);
    file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Final root file for CEA");
    if (!file_final || file_final->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    tree_final = new TTree("Events", "Final CEA tree");

    // Variables to create the new tree + linking them to the tree
    Double_t time_final;
    Double_t position_final[3];
    Double_t energy_final;
    uint8_t coinc_ijclab1;
    uint8_t coinc_ijclab2;
    uint8_t coinc_ucda;
    uint8_t coinc_ucdb;
    uint8_t coinc_ucdc;
    uint8_t coinc_ucdd;
    uint8_t chain_len;

    tree_final->Branch("time", &time_final, "time/D");
    tree_final->Branch("position", &position_final, "position[3]/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("coinc_ijclab1", &coinc_ijclab1, "coinc_ijclab1/b");
    tree_final->Branch("coinc_ijclab2", &coinc_ijclab2, "coinc_ijclab2/b");
    tree_final->Branch("coinc_ucda", &coinc_ucda, "coinc_ucda/b");
    tree_final->Branch("coinc_ucdb", &coinc_ucdb, "coinc_ucdb/b");
    tree_final->Branch("coinc_ucdc", &coinc_ucdc, "coinc_ucdc/b");
    tree_final->Branch("coinc_ucdd", &coinc_ucdd, "coinc_ucdd/b");
    tree_final->Branch("chain_len", &chain_len, "chain_len/b");

    uint64_t nentries = tree_init->GetEntries();
    // Using vectors to easily fill the tree
    std::vector<uint8_t> vec_coinc_ijclab1(nentries, 0);
    std::vector<uint8_t> vec_coinc_ijclab2(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucda(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucdb(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucdc(nentries, 0);
    std::vector<uint8_t> vec_coinc_ucdd(nentries, 0);

    for (uint64_t i = 0; i<coinc_idx_ijclab1.size(); i++) {
        vec_coinc_ijclab1[coinc_idx_ijclab1[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ijclab2.size(); i++) {
        vec_coinc_ijclab2[coinc_idx_ijclab2[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucda.size(); i++) {
        vec_coinc_ucda[coinc_idx_ucda[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucdb.size(); i++) {
        vec_coinc_ucdb[coinc_idx_ucdb[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucdc.size(); i++) {
        vec_coinc_ucdc[coinc_idx_ucdc[i]] = 1;
    }
    for (uint64_t i = 0; i<coinc_idx_ucdd.size(); i++) {
        vec_coinc_ucdd[coinc_idx_ucdd[i]] = 1;
    }

    for (uint64_t i = 0; i<nentries; i++) {
        tree_init->GetEntry(i);

        time_final = time_init;
        position_final[0] = position_init[0];
        position_final[1] = position_init[1];
        position_final[2] = position_init[2];
        energy_final = energy_init;

        coinc_ijclab1 = vec_coinc_ijclab1[i];
        coinc_ijclab2 = vec_coinc_ijclab2[i];
        coinc_ucda = vec_coinc_ucda[i];
        coinc_ucdb = vec_coinc_ucdb[i];
        coinc_ucdc = vec_coinc_ucdc[i];
        coinc_ucdd = vec_coinc_ucdd[i];
        chain_len = vec_coinc_ijclab1[i] + vec_coinc_ijclab2[i] + vec_coinc_ucda[i] + vec_coinc_ucdb[i] + vec_coinc_ucdc[i] + vec_coinc_ucdd[i];
        tree_final->Fill();
    }

    file_init->Close();
    file_final->Write();
    file_final->Close();
}

void CombineUCDSubDets(const std::string &tree_name_a, const std::string &tree_name_b, const std::string &tree_name_c, const std::string &tree_name_d, const std::string &final_tree_name) {
    // Extraction of the root tree
    file_exists(tree_name_a);
    auto file_a = new TFile(tree_name_a.c_str());
    if (!file_a || file_a->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_a = (TTree*) file_a->Get("Events");

    // Variables to read the initial tree
    int32_t pps_cpt_init_noabs_a;
    int32_t pps_cpt_corr_noabs_a;
    int32_t gps_corr_noabs_a;
    int32_t pps_cpt_corr_abs_a;
    int32_t pps_cpt_nocorr_abs_a;
    int32_t gps_corr_abs_a;
    Double_t time_corr_abs_a;
    Double_t energy_final_a;
    Double_t position_final_a[3];
    tree_a->SetBranchAddress("pps_cpt_init", &pps_cpt_init_noabs_a);
    tree_a->SetBranchAddress("pps_cpt_corr", &pps_cpt_corr_noabs_a);
    tree_a->SetBranchAddress("gps_corr", &gps_corr_noabs_a);
    tree_a->SetBranchAddress("pps_cpt_corr_abs", &pps_cpt_corr_abs_a);
    tree_a->SetBranchAddress("pps_cpt_init_abs", &pps_cpt_nocorr_abs_a);
    tree_a->SetBranchAddress("gps_corr_abs", &gps_corr_abs_a);
    tree_a->SetBranchAddress("time_corr_abs", &time_corr_abs_a);
    tree_a->SetBranchAddress("energy", &energy_final_a);
    tree_a->Branch("position", &position_final_a, "position[3]/D");

    file_exists(tree_name_b);
    auto file_b = new TFile(tree_name_b.c_str());
    if (!file_b || file_b->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_b = (TTree*) file_b->Get("Events");

    // Variables to read the initial tree
    int32_t pps_cpt_init_noabs_b;
    int32_t pps_cpt_corr_noabs_b;
    int32_t gps_corr_noabs_b;
    int32_t pps_cpt_corr_abs_b;
    int32_t pps_cpt_nocorr_abs_b;
    int32_t gps_corr_abs_b;
    Double_t time_corr_abs_b;
    Double_t energy_final_b;
    Double_t position_final_b[3];

    tree_b->SetBranchAddress("pps_cpt_init", &pps_cpt_init_noabs_b);
    tree_b->SetBranchAddress("pps_cpt_corr", &pps_cpt_corr_noabs_b);
    tree_b->SetBranchAddress("gps_corr", &gps_corr_noabs_b);
    tree_b->SetBranchAddress("pps_cpt_corr_abs", &pps_cpt_corr_abs_b);
    tree_b->SetBranchAddress("pps_cpt_init_abs", &pps_cpt_nocorr_abs_b);
    tree_b->SetBranchAddress("gps_corr_abs", &gps_corr_abs_b);
    tree_b->SetBranchAddress("time_corr_abs", &time_corr_abs_b);
    tree_b->SetBranchAddress("energy", &energy_final_b);
    tree_b->Branch("position", &position_final_b, "position[3]/D");

    file_exists(tree_name_c);
    auto file_c = new TFile(tree_name_c.c_str());
    if (!file_c || file_c->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_c = (TTree*) file_c->Get("Events");

    // Variables to read the initial tree
    int32_t pps_cpt_init_noabs_c;
    int32_t pps_cpt_corr_noabs_c;
    int32_t gps_corr_noabs_c;
    int32_t pps_cpt_corr_abs_c;
    int32_t pps_cpt_nocorr_abs_c;
    int32_t gps_corr_abs_c;
    Double_t time_corr_abs_c;
    Double_t energy_final_c;
    Double_t position_final_c[3];

    tree_c->SetBranchAddress("pps_cpt_init", &pps_cpt_init_noabs_c);
    tree_c->SetBranchAddress("pps_cpt_corr", &pps_cpt_corr_noabs_c);
    tree_c->SetBranchAddress("gps_corr", &gps_corr_noabs_c);
    tree_c->SetBranchAddress("pps_cpt_corr_abs", &pps_cpt_corr_abs_c);
    tree_c->SetBranchAddress("pps_cpt_init_abs", &pps_cpt_nocorr_abs_c);
    tree_c->SetBranchAddress("gps_corr_abs", &gps_corr_abs_c);
    tree_c->SetBranchAddress("time_corr_abs", &time_corr_abs_c);
    tree_c->SetBranchAddress("energy", &energy_final_c);
    tree_c->Branch("position", &position_final_c, "position[3]/D");

    file_exists(tree_name_d);
    auto file_d = new TFile(tree_name_d.c_str());
    if (!file_d || file_d->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_d = (TTree*) file_d->Get("Events");

    // Variables to read the initial tree
    int32_t pps_cpt_init_noabs_d;
    int32_t pps_cpt_corr_noabs_d;
    int32_t gps_corr_noabs_d;
    int32_t pps_cpt_corr_abs_d;
    int32_t pps_cpt_nocorr_abs_d;
    int32_t gps_corr_abs_d;
    Double_t time_corr_abs_d;
    Double_t energy_final_d;
    Double_t position_final_d[3];

    tree_d->SetBranchAddress("pps_cpt_init", &pps_cpt_init_noabs_d);
    tree_d->SetBranchAddress("pps_cpt_corr", &pps_cpt_corr_noabs_d);
    tree_d->SetBranchAddress("gps_corr", &gps_corr_noabs_d);
    tree_d->SetBranchAddress("pps_cpt_corr_abs", &pps_cpt_corr_abs_d);
    tree_d->SetBranchAddress("pps_cpt_init_abs", &pps_cpt_nocorr_abs_d);
    tree_d->SetBranchAddress("gps_corr_abs", &gps_corr_abs_d);
    tree_d->SetBranchAddress("time_corr_abs", &time_corr_abs_d);
    tree_d->SetBranchAddress("energy", &energy_final_d);
    tree_d->Branch("position", &position_final_d, "position[3]/D");

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    file_final = new TFile(final_tree_name.c_str(), "RECREATE", "Abs time root file for UCD");
    if (!file_final || file_final->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    tree_final = new TTree("Events", "Absolute time UCD tree");

    // Variables to create the new tree + linking them to the tree
    int32_t pps_cpt_init_noabs;
    int32_t pps_cpt_corr_noabs;
    int32_t gps_corr_noabs;
    int32_t pps_cpt_corr_abs;
    int32_t pps_cpt_nocorr_abs;
    int32_t gps_corr_abs;
    Double_t time_corr_abs;
    Double_t energy_final;
    Double_t position_final[3];
    int8_t det_id;
    tree_final->Branch("pps_cpt_init", &pps_cpt_init_noabs, "pps_cpt_init/I");
    tree_final->Branch("pps_cpt_corr", &pps_cpt_corr_noabs, "pps_cpt_corr/I");
    tree_final->Branch("gps_corr", &gps_corr_noabs, "gps_corr/I");
    tree_final->Branch("pps_cpt_corr_abs", &pps_cpt_corr_abs, "pps_cpt_corr_abs/I");
    tree_final->Branch("pps_cpt_init_abs", &pps_cpt_nocorr_abs, "pps_cpt_init_abs/I");
    tree_final->Branch("gps_corr_abs", &gps_corr_abs, "gps_corr_abs/I");
    tree_final->Branch("time_corr_abs", &time_corr_abs, "time_corr_abs/D");
    tree_final->Branch("energy", &energy_final, "energy/D");
    tree_final->Branch("position", &position_final, "position[3]/D");
    tree_final->Branch("det_id", &det_id, "det_id/b");

    // Combining the UCD trees
    uint64_t nentries_a = tree_a->GetEntries();
    uint64_t nentries_b = tree_b->GetEntries();
    uint64_t nentries_c = tree_c->GetEntries();
    uint64_t nentries_d = tree_d->GetEntries();
    std::vector<Double_t> timevec;
    timevec.reserve(nentries_a + nentries_b + nentries_c + nentries_d);

    std::cout << "UCDA number of event : " << nentries_a << std::endl;
    for (uint64_t i = 0; i<nentries_a; i++) {
        tree_a->GetEntry(i);
        timevec.push_back(time_corr_abs_a);
    }

    std::cout << "UCDB number of event : " << nentries_b << std::endl;
    for (uint64_t i = 0; i<nentries_b; i++) {
        tree_b->GetEntry(i);
        timevec.push_back(time_corr_abs_b);
    }

    std::cout << "UCDC number of event : " << nentries_c << std::endl;
    for (uint64_t i = 0; i<nentries_c; i++) {
        tree_c->GetEntry(i);
        timevec.push_back(time_corr_abs_c);
    }

    std::cout << "UCDD number of event : " << nentries_d << std::endl;
    for (uint64_t i = 0; i<nentries_d; i++) {
        tree_d->GetEntry(i);
        timevec.push_back(time_corr_abs_d);
    }

    uint64_t nentries_final = timevec.size();
    std::vector<size_t> sorted_index = sort_index(timevec);
//    uint64_t ite_count = 0;
    for (size_t sorted_ite : sorted_index){
        // Getting the entry in sorted order
        if (sorted_ite < nentries_a) {
            tree_a->GetEntry(sorted_ite);
            pps_cpt_init_noabs = pps_cpt_init_noabs_a;
            pps_cpt_corr_noabs = pps_cpt_corr_noabs_a;
            gps_corr_noabs = gps_corr_noabs_a;
            pps_cpt_corr_abs = pps_cpt_corr_abs_a;
            pps_cpt_nocorr_abs = pps_cpt_nocorr_abs_a;
            gps_corr_abs = gps_corr_abs_a;
            time_corr_abs = time_corr_abs_a;
            energy_final = energy_final_a;
            position_final[0] = position_final_a[0];
            position_final[1] = position_final_a[1];
            position_final[2] = position_final_a[2];
            det_id = 3;
        } else if (sorted_ite < nentries_a + nentries_b) {
            tree_b->GetEntry(sorted_ite - nentries_a);
            pps_cpt_init_noabs = pps_cpt_init_noabs_b;
            pps_cpt_corr_noabs = pps_cpt_corr_noabs_b;
            gps_corr_noabs = gps_corr_noabs_b;
            pps_cpt_corr_abs = pps_cpt_corr_abs_b;
            pps_cpt_nocorr_abs = pps_cpt_nocorr_abs_b;
            gps_corr_abs = gps_corr_abs_b;
            time_corr_abs = time_corr_abs_b;
            energy_final = energy_final_b;
            position_final[0] = position_final_b[0];
            position_final[1] = position_final_b[1];
            position_final[2] = position_final_b[2];
            det_id = 4;
        } else if (sorted_ite < nentries_a + nentries_b + nentries_c) {
            tree_c->GetEntry(sorted_ite - nentries_a - nentries_b);
            pps_cpt_init_noabs = pps_cpt_init_noabs_c;
            pps_cpt_corr_noabs = pps_cpt_corr_noabs_c;
            gps_corr_noabs = gps_corr_noabs_c;
            pps_cpt_corr_abs = pps_cpt_corr_abs_c;
            pps_cpt_nocorr_abs = pps_cpt_nocorr_abs_c;
            gps_corr_abs = gps_corr_abs_c;
            time_corr_abs = time_corr_abs_c;
            energy_final = energy_final_c;
            position_final[0] = position_final_c[0];
            position_final[1] = position_final_c[1];
            position_final[2] = position_final_c[2];
            det_id = 5;
        } else {
            tree_d->GetEntry(sorted_ite - nentries_a - nentries_b - nentries_c);
            pps_cpt_init_noabs = pps_cpt_init_noabs_d;
            pps_cpt_corr_noabs = pps_cpt_corr_noabs_d;
            gps_corr_noabs = gps_corr_noabs_d;
            pps_cpt_corr_abs = pps_cpt_corr_abs_d;
            pps_cpt_nocorr_abs = pps_cpt_nocorr_abs_d;
            gps_corr_abs = gps_corr_abs_d;
            time_corr_abs = time_corr_abs_d;
            energy_final = energy_final_d;
            position_final[0] = position_final_d[0];
            position_final[1] = position_final_d[1];
            position_final[2] = position_final_d[2];
            det_id = 6;
        }

        tree_final->Fill();
//        if (ite_count % 100000 == 0) {
//            std::cout << "\rOrdering the events in the new tree : entry " << ite_count << " over " << nentries_final << std::flush;
//        }
//        ite_count += 1;
    }

    file_final->Write();

//    std::cout << "size time vector vs size final tree " << timevec.size() << " - " << tree_final->GetEntries() << std::endl;


    file_a->Close();
    file_b->Close();
    file_c->Close();
    file_d->Close();
    file_final->Close();
}

//struct DetectorHit {
//    uint8_t det_id;     // Detector id
//    Double_t time;      // absolute corrected time of the trigger
//    Double_t energy;    // Energy deposit
//    Double_t pos_x;     // x position
//    Double_t pos_y;     // y position
//    Double_t pos_z;     // z position
//};

struct DetectorTreeStruct {
    TTree* tree = nullptr;
    uint64_t nentry = 0;
    Double_t time = 0;
    Double_t energy = 0;
    Double_t pos[3] = {0};
    uint8_t det_id = 0;
};

void AttributeValues(Double_t &time, Double_t &energy, Double_t &posx, Double_t &posy, Double_t &posz, char mask[7], std::vector<uint8_t> temp_dets, std::vector<size_t> temp_tree_evnt_id, std::map<uint8_t, size_t> det_id_to_tree_idx, std::array<DetectorTreeStruct, 4> detectorTrees, uint8_t det_id) {
    uint64_t tree_entry;
    uint8_t tree_number;

    auto it = std::find(temp_dets.begin(), temp_dets.end(), det_id);

    if (it != temp_dets.end()) {
        size_t temp_lists_idx = std::distance(temp_dets.begin(), it);
        if (det_id != temp_dets[temp_lists_idx]) {
            std::cerr << "Error : the searching for id in temp_dets went wrong" << std::endl;
        }
        tree_number = det_id_to_tree_idx.at(det_id);
        tree_entry = temp_tree_evnt_id[temp_lists_idx];
        detectorTrees[tree_number].tree->GetEntry(tree_entry);
        // id found in the temp_dets vector, saving the values
//        std::cout << "IN THE FUNCTION : " << std::endl;
//        std::cout << "found det_id and det_id" << static_cast<uint32_t>(temp_dets[temp_lists_idx]) << " " << static_cast<uint32_t>(det_id) << std::endl;
//        std::cout << "tree_number : " << static_cast<uint32_t>(tree_number) << " rank : "<< tree_entry << std::endl;
//        std::cout << "rank in the temp list " << temp_lists_idx << std::endl;
//        std::cout << "time " << detectorTrees[tree_number].time << " " << std::endl;
//        std::cout << "energy " << detectorTrees[tree_number].energy << std::endl;
//        std::cout << "pos_x " << detectorTrees[tree_number].pos[0] << std::endl;
//        std::cout << "pos_y " << detectorTrees[tree_number].pos[1] << std::endl;
//        std::cout << "pos_z " << detectorTrees[tree_number].pos[2] << std::endl;

        mask[det_id - 1] = 1;
        time = detectorTrees[tree_number].time;
        energy = detectorTrees[tree_number].energy;
        posx = detectorTrees[tree_number].pos[0];
        posy = detectorTrees[tree_number].pos[1];
        posz = detectorTrees[tree_number].pos[2];
    } else {
        // id not in temp_dets vector, saving default
        mask[det_id - 1] = 0;
        time = -999;
        energy = -999;
        posx = -999;
        posy = -999;
        posz = -999;
    }

}

void MakeWholeDataTree(const std::string &fname_final_tree, const std::string &fname_cea_corr_abs, const std::string &fname_dssd_corr_abs, const std::string &fname_ucd_corr_abs, const std::string &fname_maud_corr_abs) {
    // Detectors are recognized by their id : cea 1, dssd 2, ucd 3, 4, 5, 6 and maud 7

    // Extraction of the CEA root tree
    file_exists(fname_cea_corr_abs);
    auto file_cea = new TFile(fname_cea_corr_abs.c_str());
    if (!file_cea || file_cea->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    // Extraction of the DSSD root tree
    file_exists(fname_dssd_corr_abs);
    auto file_dssd = new TFile(fname_dssd_corr_abs.c_str());
    if (!file_dssd || file_dssd->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    // Extraction of the UCD root tree
    file_exists(fname_ucd_corr_abs);
    auto file_ucd = new TFile(fname_ucd_corr_abs.c_str());
    if (!file_ucd || file_ucd->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    // Extraction of the MAUD root tree
    file_exists(fname_maud_corr_abs);
    auto file_maud = new TFile(fname_maud_corr_abs.c_str());
    if (!file_maud || file_maud->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }

    std::array<DetectorTreeStruct, 4> detectorTrees;

    // Filling the structure with the trees
    detectorTrees[0].tree = (TTree*) file_cea->Get("Events");
    detectorTrees[0].tree->SetBranchAddress("time_corr_abs", &detectorTrees[0].time);
    detectorTrees[0].tree->SetBranchAddress("energy", &detectorTrees[0].energy);
    detectorTrees[0].tree->SetBranchAddress("position", &detectorTrees[0].pos);
    detectorTrees[0].tree->SetBranchAddress("det_id", &detectorTrees[0].det_id);
    detectorTrees[0].nentry = detectorTrees[0].tree->GetEntries();

    detectorTrees[1].tree = (TTree*) file_dssd->Get("Events");
    detectorTrees[1].tree->SetBranchAddress("time_corr_abs", &detectorTrees[1].time);
    detectorTrees[1].tree->SetBranchAddress("energy", &detectorTrees[1].energy);
    detectorTrees[1].tree->SetBranchAddress("position", &detectorTrees[1].pos);
    detectorTrees[1].tree->SetBranchAddress("det_id", &detectorTrees[1].det_id);
    detectorTrees[1].nentry = detectorTrees[1].tree->GetEntries();

    detectorTrees[2].tree = (TTree*) file_ucd->Get("Events");
    detectorTrees[2].tree->SetBranchAddress("time_corr_abs", &detectorTrees[2].time);
    detectorTrees[2].tree->SetBranchAddress("energy", &detectorTrees[2].energy);
    detectorTrees[2].tree->SetBranchAddress("position", &detectorTrees[2].pos);
    detectorTrees[2].tree->SetBranchAddress("det_id", &detectorTrees[2].det_id);
    detectorTrees[2].nentry = detectorTrees[2].tree->GetEntries();

    detectorTrees[3].tree = (TTree*) file_maud->Get("Events");
    detectorTrees[3].tree->SetBranchAddress("time_corr_abs", &detectorTrees[3].time);
    detectorTrees[3].tree->SetBranchAddress("energy", &detectorTrees[3].energy);
    detectorTrees[3].tree->SetBranchAddress("position", &detectorTrees[3].pos);
    detectorTrees[3].tree->SetBranchAddress("det_id", &detectorTrees[3].det_id);
    detectorTrees[3].nentry = detectorTrees[3].tree->GetEntries();

//    auto tree_cea = (TTree*) file_cea->Get("Events");

    // Declaring file_final and tree_final
    TFile* file_final = nullptr;
    TTree* tree_final = nullptr;

    file_final = new TFile(fname_final_tree.c_str(), "RECREATE", "Full data free file");
    if (!file_final || file_final->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    tree_final = new TTree("Events", "Full data tree file");

    uint8_t event_len;
    char mask[7];

//    uint8_t cea_det_id;
    Double_t cea_time;
    Double_t cea_energy;
    Double_t cea_posx;
    Double_t cea_posy;
    Double_t cea_posz;

//    uint8_t dssd_det_id;
    Double_t dssd_time;
    Double_t dssd_energy;
    Double_t dssd_posx;
    Double_t dssd_posy;
    Double_t dssd_posz;

//    uint8_t ucda_det_id;
    Double_t ucda_time;
    Double_t ucda_energy;
    Double_t ucda_posx;
    Double_t ucda_posy;
    Double_t ucda_posz;

//    uint8_t ucdb_det_id;
    Double_t ucdb_time;
    Double_t ucdb_energy;
    Double_t ucdb_posx;
    Double_t ucdb_posy;
    Double_t ucdb_posz;

//    uint8_t ucdc_det_id;
    Double_t ucdc_time;
    Double_t ucdc_energy;
    Double_t ucdc_posx;
    Double_t ucdc_posy;
    Double_t ucdc_posz;

//    uint8_t ucdd_det_id;
    Double_t ucdd_time;
    Double_t ucdd_energy;
    Double_t ucdd_posx;
    Double_t ucdd_posy;
    Double_t ucdd_posz;

//    uint8_t maud_det_id;
    Double_t maud_time;
    Double_t maud_energy;
    Double_t maud_posx;
    Double_t maud_posy;
    Double_t maud_posz;

//    std::array<Double_t, 5> cea_vals;
//    std::array<Double_t, 5> dssd_vals;
//    std::array<Double_t, 5> ucda_vals;
//    std::array<Double_t, 5> ucdb_vals;
//    std::array<Double_t, 5> ucdc_vals;
//    std::array<Double_t, 5> ucdd_vals;
//    std::array<Double_t, 5> maud_vals;

    tree_final->Branch("event_len", &event_len, "event_len/b");
//    tree_final->Branch("mask", &mask, "mask[7]/b");

////    tree_final->Branch("cea_det_id", &cea_det_id, "cea_det_id/b");
//    tree_final->Branch("cea_time", &cea_vals[0], "cea_time/D");
//    tree_final->Branch("cea_energy", &cea_vals[1], "cea_energy/D");
//    tree_final->Branch("cea_posx", &cea_vals[2], "cea_posx/D");
//    tree_final->Branch("cea_posy", &cea_vals[3], "cea_posy/D");
//    tree_final->Branch("cea_posz", &cea_vals[4], "cea_posz/D");
//
////    tree_final->Branch("dssd_det_id", &dssd_det_id, "dssd_det_id/b");
//    tree_final->Branch("dssd_time", &dssd_vals[0], "dssd_time/D");
//    tree_final->Branch("dssd_energy", &dssd_vals[1], "dssd_energy/D");
//    tree_final->Branch("dssd_posx", &dssd_vals[2], "dssd_posx/D");
//    tree_final->Branch("dssd_posy", &dssd_vals[3], "dssd_posy/D");
//    tree_final->Branch("dssd_posz", &dssd_vals[4], "dssd_posz/D");
//
////    tree_final->Branch("ucda_det_id", &ucda_det_id, "ucda_det_id/b");
//    tree_final->Branch("ucda_time", &ucda_vals[0], "ucda_time/D");
//    tree_final->Branch("ucda_energy", &ucda_vals[1], "ucda_energy/D");
//    tree_final->Branch("ucda_posx", &ucda_vals[2], "ucda_posx/D");
//    tree_final->Branch("ucda_posy", &ucda_vals[3], "ucda_posy/D");
//    tree_final->Branch("ucda_posz", &ucda_vals[4], "ucda_posz/D");
//
////    tree_final->Branch("ucdb_det_id", &ucdb_det_id, "ucdb_det_id/b");
//    tree_final->Branch("ucdb_time", &ucdb_vals[0], "ucdb_time/D");
//    tree_final->Branch("ucdb_energy", &ucdb_vals[1], "ucdb_energy/D");
//    tree_final->Branch("ucdb_posx", &ucdb_vals[2], "ucdb_posx/D");
//    tree_final->Branch("ucdb_posy", &ucdb_vals[3], "ucdb_posy/D");
//    tree_final->Branch("ucdb_posz", &ucdb_vals[4], "ucdb_posz/D");
//
////    tree_final->Branch("ucdc_det_id", &ucdc_det_id, "ucdc_det_id/b");
//    tree_final->Branch("ucdc_time", &ucdc_vals[0], "ucdc_time/D");
//    tree_final->Branch("ucdc_energy", &ucdc_vals[1], "ucdc_energy/D");
//    tree_final->Branch("ucdc_posx", &ucdc_vals[2], "ucdc_posx/D");
//    tree_final->Branch("ucdc_posy", &ucdc_vals[3], "ucdc_posy/D");
//    tree_final->Branch("ucdc_posz", &ucdc_vals[4], "ucdc_posz/D");
//
////    tree_final->Branch("ucdd_det_id", &ucdd_det_id, "ucdd_det_id/b");
//    tree_final->Branch("ucdd_time", &ucdd_vals[0], "ucdd_time/D");
//    tree_final->Branch("ucdd_energy", &ucdd_vals[1], "ucdd_energy/D");
//    tree_final->Branch("ucdd_posx", &ucdd_vals[2], "ucdd_posx/D");
//    tree_final->Branch("ucdd_posy", &ucdd_vals[3], "ucdd_posy/D");
//    tree_final->Branch("ucdd_posz", &ucdd_vals[4], "ucdd_posz/D");

////    tree_final->Branch("maud_det_id", &maud_det_id, "maud_det_id/b");
//    tree_final->Branch("maud_time", &maud_vals[0], "maud_time/D");
//    tree_final->Branch("maud_energy", &maud_vals[1], "maud_energy/D");
//    tree_final->Branch("maud_posx", &maud_vals[2], "maud_posx/D");
//    tree_final->Branch("maud_posy", &maud_vals[3], "maud_posy/D");
//    tree_final->Branch("maud_posz", &maud_vals[4], "maud_posz/D");

//    tree_final->Branch("cea_det_id", &cea_det_id, "cea_det_id/b");
    tree_final->Branch("cea_time", &cea_time, "cea_time/D");
    tree_final->Branch("cea_energy", &cea_energy, "cea_energy/D");
    tree_final->Branch("cea_posx", &cea_posx, "cea_posx/D");
    tree_final->Branch("cea_posy", &cea_posy, "cea_posy/D");
    tree_final->Branch("cea_posz", &cea_posz, "cea_posz/D");

//    tree_final->Branch("dssd_det_id", &dssd_det_id, "dssd_det_id/b");
    tree_final->Branch("dssd_time", &dssd_time, "dssd_time/D");
    tree_final->Branch("dssd_energy", &dssd_energy, "dssd_energy/D");
    tree_final->Branch("dssd_posx", &dssd_posx, "dssd_posx/D");
    tree_final->Branch("dssd_posy", &dssd_posy, "dssd_posy/D");
    tree_final->Branch("dssd_posz", &dssd_posz, "dssd_posz/D");

//    tree_final->Branch("ucda_det_id", &ucda_det_id, "ucda_det_id/b");
    tree_final->Branch("ucda_time", &ucda_time, "ucda_time/D");
    tree_final->Branch("ucda_energy", &ucda_energy, "ucda_energy/D");
    tree_final->Branch("ucda_posx", &ucda_posx, "ucda_posx/D");
    tree_final->Branch("ucda_posy", &ucda_posy, "ucda_posy/D");
    tree_final->Branch("ucda_posz", &ucda_posz, "ucda_posz/D");

//    tree_final->Branch("ucdb_det_id", &ucdb_det_id, "ucdb_det_id/b");
    tree_final->Branch("ucdb_time", &ucdb_time, "ucdb_time/D");
    tree_final->Branch("ucdb_energy", &ucdb_energy, "ucdb_energy/D");
    tree_final->Branch("ucdb_posx", &ucdb_posx, "ucdb_posx/D");
    tree_final->Branch("ucdb_posy", &ucdb_posy, "ucdb_posy/D");
    tree_final->Branch("ucdb_posz", &ucdb_posz, "ucdb_posz/D");

//    tree_final->Branch("ucdc_det_id", &ucdc_det_id, "ucdc_det_id/b");
    tree_final->Branch("ucdc_time", &ucdc_time, "ucdc_time/D");
    tree_final->Branch("ucdc_energy", &ucdc_energy, "ucdc_energy/D");
    tree_final->Branch("ucdc_posx", &ucdc_posx, "ucdc_posx/D");
    tree_final->Branch("ucdc_posy", &ucdc_posy, "ucdc_posy/D");
    tree_final->Branch("ucdc_posz", &ucdc_posz, "ucdc_posz/D");

//    tree_final->Branch("ucdd_det_id", &ucdd_det_id, "ucdd_det_id/b");
    tree_final->Branch("ucdd_time", &ucdd_time, "ucdd_time/D");
    tree_final->Branch("ucdd_energy", &ucdd_energy, "ucdd_energy/D");
    tree_final->Branch("ucdd_posx", &ucdd_posx, "ucdd_posx/D");
    tree_final->Branch("ucdd_posy", &ucdd_posy, "ucdd_posy/D");
    tree_final->Branch("ucdd_posz", &ucdd_posz, "ucdd_posz/D");

//    tree_final->Branch("maud_det_id", &maud_det_id, "maud_det_id/b");
    tree_final->Branch("maud_time", &maud_time, "maud_time/D");
    tree_final->Branch("maud_energy", &maud_energy, "maud_energy/D");
    tree_final->Branch("maud_posx", &maud_posx, "maud_posx/D");
    tree_final->Branch("maud_posy", &maud_posy, "maud_posy/D");
    tree_final->Branch("maud_posz", &maud_posz, "maud_posz/D");

//    std::array<std::reference_wrapper<std::array<Double_t, 5>>, 7> event_values = {
//    cea_vals, dssd_vals, ucda_vals, ucdb_vals, ucdc_vals, ucdd_vals, maud_vals};

    // Combining the UCD trees
    std::vector<Double_t> timevec;
    std::vector<uint64_t> det_id;
    std::vector<uint64_t> det_tree_index;
    timevec.reserve(detectorTrees[0].nentry + detectorTrees[1].nentry + detectorTrees[2].nentry + detectorTrees[3].nentry);

    for (size_t ite_det = 0; ite_det<4; ++ite_det) {
//        std::cout << "DET : " << ite_det << std::endl;
        for (uint64_t i = 0; i<detectorTrees[ite_det].nentry; ++i) {
            detectorTrees[ite_det].tree->GetEntry(i);
            timevec.push_back(detectorTrees[ite_det].time);
            det_id.push_back(detectorTrees[ite_det].det_id);
            det_tree_index.push_back(i);
//            if ( i == 0 || i == 1 || i == 2 || i == 3) {
////                OK
//                std::cout << "Show first ev in dets" << std::endl;
//                std::cout << "time " << detectorTrees[ite_det].time << std::endl;
//                std::cout << "energy " << detectorTrees[ite_det].energy << std::endl;
//                std::cout << "x " << detectorTrees[ite_det].pos[0] << std::endl;
//                std::cout << "y " << detectorTrees[ite_det].pos[1] << std::endl;
//                std::cout << "z " << detectorTrees[ite_det].pos[2] << std::endl;
//                std::cout << "det_id, " << detectorTrees[ite_det].det_id << std::endl;
//                std::cout << "ite det, " << ite_det << std::endl;
//                std::cout << "det_tree_index, " << i << std::endl;
//            }
        }
    }

    uint64_t nentries_temp = timevec.size();
    std::vector<size_t> sorted_index = sort_index(timevec);

    uint64_t i_loop = 0;
    uint64_t j_loop = 0;
    uint64_t idx_i;
    uint64_t idx_j;
    Double_t window = 2e-6;
    std::vector<uint8_t> temp_dets;
    std::vector<size_t> temp_tree_evnt_id;
    std::unordered_set<size_t> seen;

    std::map<uint8_t, size_t> det_id_to_tree_idx = {
    {1, 0}, {2, 1}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 3}};

    uint64_t coincwarning = 0;

    // Reading first entry of trees
    // Necessary for a proper functioning, but no idea why
    detectorTrees[0].tree->GetEntry(0);
    detectorTrees[1].tree->GetEntry(0);
    detectorTrees[2].tree->GetEntry(0);
    detectorTrees[3].tree->GetEntry(0);

//    i_loop = 2000000;
//    j_loop = 2000000;
//    nentries_temp = 3000000;
    while (i_loop < nentries_temp) {
//        std::cout << "\nNew loop" << std::endl;
        temp_dets.clear();
        temp_tree_evnt_id.clear();
        seen.clear();
        idx_i = sorted_index[i_loop];
        idx_j = sorted_index[j_loop];
//        std::cout << "i et j: " << i_loop << j_loop << std::endl;
//        std::cout << "itri et jtri: " << idx_i << idx_j << std::endl;
//        std::cout << "censés être les mêmes" << std::endl;
        if (timevec[idx_j] < timevec[idx_i]) {
            std::cerr << "WARNING, the time vector is supposed to be ordered but i+1 val < i val" <<std::endl;
        }
//        std::cout << "BEGIN : "<< std::endl;
        while (j_loop < nentries_temp && timevec[idx_j] - timevec[idx_i] <= window) {
//            std::cout << "  time values : " << timevec[idx_i] << "  " << timevec[idx_j] << std::endl;
//            if (timevec[idx_j] - timevec[idx_i] == 0) {
//                std::cout << " SAME  diff time values : " << timevec[idx_j] - timevec[idx_i] << std::endl;
//            } else {
//                std::cout << " NOT SAME  diff time values : " << timevec[idx_j] - timevec[idx_i] << std::endl;
//            }
            // Handling the case where the same detector triggered twice in the same time window.
            // Probably to optimize, here we just suppress one of the informations maybe the wrong one is deleted
            if (seen.insert(det_id[idx_j]).second) {
//                std::cout << "adding unique det" << det_id[idx_j] << std::endl;
//                std::cout << "adding unique idx" << idx_j << std::endl;
                // Stores the id of the detector
                temp_dets.push_back(det_id[idx_j]);
                // Stores the index of the event in the detectors tree
                temp_tree_evnt_id.push_back(det_tree_index[idx_j]);
//                std::cout << " === i and j  : " << i_loop << " " << j_loop << std::endl;
//                std::cout << " === idi and idj  : " << idx_i << " " << idx_j << std::endl;
//                std::cout << " === times i and j  : " << timevec[idx_i] << " " << timevec[idx_j] << std::endl;
//                std::cout << " === det i and j : " << det_id[idx_i]  << " " << det_id[idx_j]<< std::endl;
//                std::cout << " === tree number i and j : " << det_id_to_tree_idx[det_id[idx_i]] << " " << det_id_to_tree_idx[det_id[idx_j]] << std::endl;
//                std::cout << " === number of tree event i and j : " << det_tree_index[idx_i] << " " << det_tree_index[idx_j] << std::endl;
                detectorTrees[det_id_to_tree_idx[det_id[idx_j]]].tree->GetEntry(det_tree_index[idx_j]);
//                if (det_id[idx_j] == 3) {
//                    std::cout << "ucda tree time : " << detectorTrees[det_id_to_tree_idx[det_id[idx_j]]].time << std::endl;
//                } else if (det_id[idx_j] == 4) {
//                    std::cout << "ucdb tree time : " << detectorTrees[det_id_to_tree_idx[det_id[idx_j]]].time << std::endl;
//                } else if (det_id[idx_j] == 5) {
//                    std::cout << "ucdc tree time : " << detectorTrees[det_id_to_tree_idx[det_id[idx_j]]].time << std::endl;
//                } else if (det_id[idx_j] == 6) {
//                    std::cout << "ucdd tree time : " << detectorTrees[det_id_to_tree_idx[det_id[idx_j]]].time << std::endl;
//                }
            } else {
                coincwarning += 1;
            }

            j_loop++;
            if (j_loop < nentries_temp) {
                idx_j = sorted_index[j_loop];
            }
        }

//        std::cout << "Looping aver coincident events done : i and j : " << i_loop << j_loop << std::endl;
        i_loop = j_loop;
        // Filling the root tree
        event_len = temp_tree_evnt_id.size();

//        std::cout << std::fixed << std::setprecision(10);
        // cea vals
        AttributeValues(cea_time, cea_energy, cea_posx, cea_posy, cea_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 1);
        // dssd vals
        AttributeValues(dssd_time, dssd_energy, dssd_posx, dssd_posy, dssd_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 2);
        // ucda vals
        AttributeValues(ucda_time, ucda_energy, ucda_posx, ucda_posy, ucda_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 3);
        // ucdb vals
        AttributeValues(ucdb_time, ucdb_energy, ucdb_posx, ucdb_posy, ucdb_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 4);
        // ucdc vals
        AttributeValues(ucdc_time, ucdc_energy, ucdc_posx, ucdc_posy, ucdc_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 5);
        // ucdd vals
        AttributeValues(ucdd_time, ucdd_energy, ucdd_posx, ucdd_posy, ucdd_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 6);
        // maud vals
        AttributeValues(maud_time, maud_energy, maud_posx, maud_posy, maud_posz, mask, temp_dets, temp_tree_evnt_id, det_id_to_tree_idx, detectorTrees, 7);
//        std::cout << "len : " << temp_dets.size() << std::endl;
//        std::cout << "times : " << cea_time << " " << dssd_time << " " << ucda_time << " " << ucdb_time << " " << ucdc_time << " " << ucdd_time << " " << maud_time << " " << std::endl;
//        std::cout << "energies : " << cea_energy << " " << dssd_energy << " " << ucda_energy << " " << ucdb_energy << " " << ucdc_energy << " " << ucdd_energy << " " << maud_energy << " " << std::endl;
//        std::cout << "\n" << std::endl;
        tree_final->Fill();
    }
    file_final->Write();
    file_final->Close();
    file_cea->Close();
    file_dssd->Close();
    file_ucd->Close();
    file_maud->Close();
    if (coincwarning > 0) {
        std::cout << "WARNING : There is " << coincwarning << "coincidences between events inside the same detector" << std::endl;
    }
}

void FindCoincidences(const std::string &tree_name_a, const std::string &tree_name_b, std::vector<uint64_t> &idx_tree_a, std::vector<uint64_t> &idx_tree_b, std::vector<Double_t> &coinc_delays, const Double_t window_delay, const bool display, std::vector<TCanvas*> &save_canvas, const std::string &comment) {
    // Open the 2 trees and link the branches
    auto file_a = new TFile(tree_name_a.c_str());
    if (!file_a || file_a->IsZombie()) {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
    auto tree_a = (TTree*) file_a->Get("Events");

    auto file_b = new TFile(tree_name_b.c_str());
    if (!file_b || file_b->IsZombie()) {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
    auto tree_b = (TTree*) file_b->Get("Events");

    Double_t time_a;
    tree_a->SetBranchAddress("time_corr_abs",     &time_a);
    Double_t time_b;
    tree_b->SetBranchAddress("time_corr_abs",     &time_b);

    // Building 2 vectors to hold the times, not to use the full tree for the operations
    uint64_t nentry_a = tree_a->GetEntries();
    uint64_t nentry_b = tree_b->GetEntries();
    std::vector<Double_t> times_a(nentry_a);
    std::vector<Double_t> times_b(nentry_b);
    for (uint64_t i = 0; i<nentry_a; i++) {
        tree_a->GetEntry(i);
        times_a[i] = time_a;
    }
    for (uint64_t i = 0; i<nentry_b; i++) {
        tree_b->GetEntry(i);
        times_b[i] = time_b;
    }

    uint64_t coinc_counter_delay = 0;
    uint64_t j_restart_delay = 0;
    bool first_ev_in_win_delay;

    for (uint64_t i = 0; i<nentry_a; i++) {
        first_ev_in_win_delay = true;
        for (uint64_t j = j_restart_delay; j<nentry_b; j++) {
            if (InWindow(times_a[i], window_delay / 2, times_b[j])) {
                coinc_counter_delay += 1;
                coinc_delays.push_back(times_a[i] - times_b[j]);
                idx_tree_a.push_back(i);
                idx_tree_b.push_back(j);
                if (first_ev_in_win_delay) {
                    first_ev_in_win_delay = false;
                    j_restart_delay = j;
                }
            }
            if (times_b[j] > times_a[i] + window_delay / 2) {
                break;
            }
        }
    }
    std::cout << "  Number of coincidences for a window of " << window_delay << " s is : " << coinc_counter_delay << std::endl;
    std::cout << "    Over a total of " << nentry_a << " events (" << coinc_counter_delay / nentry_a * 100 << "%) - " << tree_name_a << std::endl;
    std::cout << "    Over a total of " << nentry_b << " events (" << coinc_counter_delay / nentry_b * 100 << "%) - " << tree_name_b << std::endl;
    if (display) {
        auto can5 = new TCanvas((std::string("can-") + comment).c_str(), (std::string("Coinc delay histogram ") + comment).c_str(), 900, 600);
        int bin_num = window_delay / 2 / 2e-8;
        std::vector<double> weights(coinc_delays.size(),1);
        std::cout << "n bins : " << bin_num << std::endl;
        std::cout << "size coinc list : " << coinc_delays.size() << std::endl;

        can5->SetLogy();
        can5->cd(1);
        auto h11 = new TH1D((std::string("hist-") + comment + std::string("-delays")).c_str(), (std::string("Coinc delay histogram ") + comment).c_str(), bin_num, -window_delay / 2,  window_delay / 2);
        h11->FillN(coinc_delays.size(), coinc_delays.data(), weights.data());
        h11->GetXaxis()->SetTitleSize(0.045);
        h11->GetYaxis()->SetTitleSize(0.045);
        h11->GetXaxis()->SetLabelSize(0.04);
        h11->GetYaxis()->SetLabelSize(0.04);
        h11->GetXaxis()->SetTitle((std::string("Delay time ") + comment + std::string(" (s)")).c_str());
        h11->GetYaxis()->SetTitle("Number of coincidences");
        h11->SetTitle("");
        h11->Draw();

        can5->Update();
        save_canvas.push_back(can5);
        can5->SaveAs((std::string("/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/coinc/hist-") + comment + std::string(".png")).c_str());
    }
}

void DisplayDSSDChannels(const std::string &tree_name) {
    // Extraction of the root tree
    file_exists(tree_name);
    auto tree_file = new TFile(tree_name.c_str());
    if (!tree_file || tree_file->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_hist = (TTree*) tree_file->Get("Events");

    // Variables to read the initial tree
    Double_t time_corr_abs;
    uint16_t adc_channels[64];
    Double_t adc_channels_corr[64];
    Double_t adc_channels_calib[64];
    Double_t temperature;
    tree_hist->SetBranchAddress("time_corr_abs",     &time_corr_abs);
    tree_hist->SetBranchAddress("adc_channels",      &adc_channels);
    tree_hist->SetBranchAddress("adc_channels_corr",      &adc_channels_corr);
    tree_hist->SetBranchAddress("adc_channels_calib",      &adc_channels_calib);
    tree_hist->SetBranchAddress("temperature", &temperature);

    uint64_t nentries = tree_hist->GetEntries();
    uint32_t init_time;
    uint32_t end_time;
    tree_hist->GetEntry(0);
    init_time = static_cast<uint32_t>(time_corr_abs);
    tree_hist->GetEntry(nentries - 1);
    end_time = static_cast<uint32_t>(time_corr_abs + 1) ;


    TApplication DisplayApp("App", nullptr, nullptr);

    int xlimmin = 0;
    int xlimmax = 1500;
    int xlimmincalib = 0;
    int xlimmaxcalib = 1500;

////  =============== 1D spectra ===============
//    std::cout << "Filling 1D spectra" << std::endl;
////    auto can_p = new TCanvas("Spectra_p", "Spectra for p channels");
//    auto can_n = new TCanvas("Spectra_n", "Spectra for n channels");
//    can_n->Divide(8, 4);
//
//    std::vector<TH1I*> spectra(32);
//    for (size_t i=0; i<32; i++) {
//        spectra[i] = new TH1I((std::string("spec") + std::to_string(i)).c_str(), "adc spectrum", 1024, 0, 1024);
//    }
//
//    for (uint64_t i = 0; i < nentries; ++i) {
//        tree_hist->GetEntry(i);
////        std::cout << "Entry : " << time_corr_abs << " OK " << std::endl;
//        for (size_t j=0; j<32; j++) {
//        spectra[j]->Fill(adc_channels[j + 32]);
//        }
//    }
//
//    for (size_t i=0; i<32; i++) {
//        can_n->cd(i + 1);
//        spectra[i]->SetStats(0);
//        spectra[i]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
//        spectra[i]->Draw();
//        gPad->SetLogy();
//    }

////  =============== 1D spectra calibrated ===============
//    auto can_n_calib = new TCanvas("Spectra_calibrated_n", "Spectra calibrated for n channels");
//    can_n_calib->Divide(8, 4);
//
//    Double_t fitxmin = 0;
//    Double_t fitxmax = 1500;
//
//    std::vector<TH1D*> spectra_calib(32);
//    std::vector<TF1*> spec_calib_fit(32);
//    for (size_t i=0; i<32; i++) {
//        spectra_calib[i] = new TH1D((std::string("spec_calib") + std::to_string(i)).c_str(), "Energy spectrum", 600, 0, 1200);
////        spec_calib_fit[i] = new TF1((std::string("gauss_fit") + std::to_string(i)).c_str(), "gaus", fitxmin, fitxmax);
////        spec_calib_fit[i]->SetLineColor(kRed);
//    }
//
//    for (uint64_t i = 0; i < nentries; ++i) {
//        tree_hist->GetEntry(i);
////        std::cout << "Entry : " << time_corr_abs << " OK " << std::endl;
//        for (size_t j=0; j<32; j++) {
//            spectra_calib[j]->Fill(adc_channels_calib[j + 32]);
//        }
//    }
//
//    for (size_t i=0; i<32; i++) {
//        can_n_calib->cd(i + 1);
//        spectra_calib[i]->SetStats(0);
//        spectra_calib[i]->GetXaxis()->SetRangeUser(xlimmincalib, xlimmaxcalib);
////        spectra_calib[i]->Fit((std::string("gauss_fit") + std::to_string(i)).c_str(), "", "", fitxmin, fitxmax);
//        spectra_calib[i]->Draw();
//        gPad->SetLogy();
//    }


//  =============== time vs temperature ===============
//    auto can_time_temp = new TCanvas("time_vs_temp", "Time vs temperature");
//    tree_hist->Draw("temperature:time_corr_abs");

    bool show_fit_bulk = false;
    bool show_fit_peak = false;
    Double_t fitxmin2 = 450;
    Double_t fitxmax2 = 700;

    std::vector<Double_t> fitminpedest = {70, 60, 60, 120, 50, 60, 120, 55, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 110, 160, 50, 50, 125, 180, 190, 65};
    std::vector<Double_t> fitmaxpedest = {105, 100, 95, 160, 100, 90, 150, 90, 100, 95, 95, 90, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 155, 290, 95, 95, 165, 320, 220, 100};
    std::vector<Double_t> fitminbulk = {480, 500, 490, 480, 490, 480, 480, 470,  480, 480, 480, 480, 480, 480, 470, 480,  480, 470, 470, 480, 480, 480, 480, 480,  470, 470, 480, 480, 480, 480, 500, 480};
    std::vector<Double_t> fitmaxbulk = {700, 720, 710, 690, 700, 690, 690, 690,  690, 690, 690, 690, 690, 690, 690, 690,  690, 690, 690, 690, 690, 690, 690, 690,  690, 690, 700, 690, 690, 690, 710, 690};
    std::vector<Double_t> fitminpedestcorr = {70, 60, 60, 120, 50, 60, 120, 55, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 110, 160, 50, 50, 125, 180, 190, 65};
    std::vector<Double_t> fitmaxpedestcorr = {105, 100, 95, 160, 100, 90, 150, 90, 100, 95, 95, 90, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 155, 290, 95, 95, 165, 320, 220, 100};
    std::vector<Double_t> fitminbulkcorr = {480, 500, 490, 480, 490, 480, 480, 470,  480, 480, 480, 480, 480, 480, 470, 480,  480, 470, 470, 480, 480, 480, 480, 480,  470, 470, 480, 480, 480, 480, 500, 480};
    std::vector<Double_t> fitmaxbulkcorr = {700, 720, 710, 690, 700, 690, 690, 690,  690, 690, 690, 690, 690, 690, 690, 690,  690, 690, 690, 690, 690, 690, 690, 690,  690, 690, 700, 690, 690, 690, 710, 690};

    // Spec histograms
    std::vector<std::vector<TH1D*>> spectra_calib_timed(6, std::vector<TH1D*>(32));
    std::vector<std::vector<TH1D*>> spectra_corr_timed(6, std::vector<TH1D*>(32));
    std::vector<std::vector<TH1D*>> spectra_raw_timed(6, std::vector<TH1D*>(32));
    // Spec fits
    std::vector<std::vector<TF1*>> spec_raw_fit(6, std::vector<TF1*>(32));
    std::vector<std::vector<TF1*>> spec_raw_fit_pedest(6, std::vector<TF1*>(32));
    std::vector<std::vector<TF1*>> spec_corr_fit(6, std::vector<TF1*>(32));
    std::vector<std::vector<TF1*>> spec_corr_fit_pedest(6, std::vector<TF1*>(32));
    // Fit raw channels
    std::vector<std::vector<Double_t>> fit_pedest(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_pedest_err(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_bulk(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_bulk_err(6, std::vector<Double_t>(32));
    // Fit raw corrected channels
    std::vector<std::vector<Double_t>> fit_corr_pedest(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_corr_pedest_err(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_corr_bulk(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_corr_bulk_err(6, std::vector<Double_t>(32));
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            // Spec histograms
            spectra_calib_timed[i][j] = new TH1D((std::string("spec_calib_timed") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "Energy spectrum", 750, 0, 1500);
            spectra_corr_timed[i][j] = new TH1D((std::string("spec_corr_timed") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "Energy corr spectrum", 750, 0, 1500);
            spectra_raw_timed[i][j] = new TH1D((std::string("spec_raw_timed") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "Raw energy spectrum", 1024, 0, 1024);
            // Initializing the raw fits
            spec_raw_fit[i][j] = new TF1((std::string("gauss_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminbulk[j], fitmaxbulk[j]);
            spec_raw_fit[i][j]->SetLineColor(kRed);
            spec_raw_fit_pedest[i][j] = new TF1((std::string("pedest_gauss_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminpedest[j], fitmaxpedest[j]);
            spec_raw_fit_pedest[i][j]->SetLineColor(kRed);
            // Initializing the corrected raw fits
            spec_corr_fit[i][j] = new TF1((std::string("gauss_corr_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminbulkcorr[j], fitmaxbulkcorr[j]);
            spec_corr_fit[i][j]->SetLineColor(kRed);
            spec_corr_fit_pedest[i][j] = new TF1((std::string("pedest_gauss_corr_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminpedestcorr[j], fitmaxpedestcorr[j]);
            spec_corr_fit_pedest[i][j]->SetLineColor(kRed);
        }
    }

    // Filling histograms
    for (uint64_t i = 0; i < nentries; ++i) {
        tree_hist->GetEntry(i);
        if (temperature >= 10) {
            for (size_t j=0; j<32; j++) {
                spectra_calib_timed[0][j]->Fill(adc_channels_calib[j + 32]);
                spectra_corr_timed[0][j]->Fill(adc_channels_corr[j + 32]);
                spectra_raw_timed[0][j]->Fill(adc_channels[j + 32]);
            }
        } else if (temperature < 10 && temperature >= 5) {
            for (size_t j=0; j<32; j++) {
                spectra_calib_timed[1][j]->Fill(adc_channels_calib[j + 32]);
                spectra_corr_timed[1][j]->Fill(adc_channels_corr[j + 32]);
                spectra_raw_timed[1][j]->Fill(adc_channels[j + 32]);
            }
        } else if (temperature < 5 && temperature >= 0) {
            for (size_t j=0; j<32; j++) {
                spectra_calib_timed[2][j]->Fill(adc_channels_calib[j + 32]);
                spectra_corr_timed[2][j]->Fill(adc_channels_corr[j + 32]);
                spectra_raw_timed[2][j]->Fill(adc_channels[j + 32]);
            }
        } else if (temperature < 0 && temperature >= -5) {
            for (size_t j=0; j<32; j++) {
                spectra_calib_timed[3][j]->Fill(adc_channels_calib[j + 32]);
                spectra_corr_timed[3][j]->Fill(adc_channels_corr[j + 32]);
                spectra_raw_timed[3][j]->Fill(adc_channels[j + 32]);
            }
        } else if (temperature < -5 && temperature >= -10) {
            for (size_t j=0; j<32; j++) {
                spectra_calib_timed[4][j]->Fill(adc_channels_calib[j + 32]);
                spectra_corr_timed[4][j]->Fill(adc_channels_corr[j + 32]);
                spectra_raw_timed[4][j]->Fill(adc_channels[j + 32]);
            }
        } else if (temperature < -10) {
            for (size_t j=0; j<32; j++) {
                spectra_calib_timed[5][j]->Fill(adc_channels_calib[j + 32]);
                spectra_corr_timed[5][j]->Fill(adc_channels_corr[j + 32]);
                spectra_raw_timed[5][j]->Fill(adc_channels[j + 32]);
            }
        }
    }

    // Colors for different energy ranges
    int colors[6] = {kRed, kBlue, kGreen+2, kMagenta, kOrange, kCyan};

    // ============================================================================
    // Calibrated spectra
    // ============================================================================
    auto can_n_calib_timed = new TCanvas("Spectra_calibrated_n_timed", "Spectra calibrated");

    // Legend area
    TPad *pad_legend_calib = new TPad("pad_legend_calib", "Legend Pad", 0, 0.92, 1, 1);
    pad_legend_calib->SetFillColor(0);
    pad_legend_calib->Draw();
    pad_legend_calib->cd();
    pad_legend_calib->Divide(3, 1);
    std::vector<TLegend*> legend_calib(3);
    pad_legend_calib->cd(1);
    legend_calib[0] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_calib[0]->AddEntry(spectra_calib_timed[0][0], "n channels >= 10#circC", "l");
    legend_calib[0]->AddEntry(spectra_calib_timed[1][0], "n channels 10 - 5#circC", "l");
    legend_calib[0]->Draw();
    pad_legend_calib->cd(2);
    legend_calib[1] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_calib[1]->AddEntry(spectra_calib_timed[2][0], "n channels 5 - 0#circC", "l");
    legend_calib[1]->AddEntry(spectra_calib_timed[3][0], "n channels 0 - -5#circC", "l");
    legend_calib[1]->Draw();
    pad_legend_calib->cd(3);
    legend_calib[2] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_calib[2]->AddEntry(spectra_calib_timed[4][0], "n channels -5 - -10#circC", "l");
    legend_calib[2]->AddEntry(spectra_calib_timed[5][0], "n channels < -10#circC", "l");
    legend_calib[2]->Draw();

    // Histograms area
    can_n_calib_timed->cd();
    TPad *pad_hist_calib = new TPad("pad_hist_calib", "Histogram pad", 0, 0, 1, 0.92);
    pad_hist_calib->Draw();
    pad_hist_calib->cd();
    pad_hist_calib->Divide(8, 4, 0.001, 0.001);

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            pad_hist_calib->cd(j + 1);
            spectra_calib_timed[i][j]->SetStats(0);
            spectra_calib_timed[i][j]->GetXaxis()->SetRangeUser(xlimmincalib, xlimmaxcalib);
            spectra_calib_timed[i][j]->SetLineColor(colors[i]);
//            spectra_calib_timed[i][j]->Fit((std::string("gauss_fit_timed1") + std::to_string(i)).c_str(), "", "", fitxmin2, fitxmax2);
            if (i == 0) {
                spectra_calib_timed[i][j]->Draw();
            } else {
                spectra_calib_timed[i][j]->DrawCopy("same");
            }
            gPad->SetLogy();
        }
    }

    can_n_calib_timed->Update();

    // ============================================================================
    // Raw spectra
    // ============================================================================
    auto can_n_raw_timed = new TCanvas("Spectra_raw_n_timed", "Raw spectra");

    // Legend area
    TPad *pad_legend_raw = new TPad("pad_legend_raw", "Legend Pad", 0, 0.92, 1, 1);
    pad_legend_raw->SetFillColor(0);
    pad_legend_raw->Draw();
    pad_legend_raw->cd();
    pad_legend_raw->Divide(3, 1);
    std::vector<TLegend*> legend_raw(3);
    pad_legend_raw->cd(1);
    legend_raw[0] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_raw[0]->AddEntry(spectra_raw_timed[0][0], "n channels >= 10#circC", "l");
    legend_raw[0]->AddEntry(spectra_raw_timed[1][0], "n channels 10 - 5#circC", "l");
    legend_raw[0]->Draw();
    pad_legend_raw->cd(2);
    legend_raw[1] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_raw[1]->AddEntry(spectra_raw_timed[2][0], "n channels 5 - 0#circC", "l");
    legend_raw[1]->AddEntry(spectra_raw_timed[3][0], "n channels 0 - -5#circC", "l");
    legend_raw[1]->Draw();
    pad_legend_raw->cd(3);
    legend_raw[2] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_raw[2]->AddEntry(spectra_raw_timed[4][0], "n channels -5 - -10#circC", "l");
    legend_raw[2]->AddEntry(spectra_raw_timed[5][0], "n channels < -10#circC", "l");
    legend_raw[2]->Draw();

    // Histograms area
    can_n_raw_timed->cd();
    TPad *pad_hist_raw = new TPad("pad_hist_raw", "Histogram pad", 0, 0, 1, 0.92);
    pad_hist_raw->Draw();
    pad_hist_raw->cd();
    pad_hist_raw->Divide(8, 4, 0.001, 0.001);

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            pad_hist_raw->cd(j + 1);
            spectra_raw_timed[i][j]->SetStats(0);
            spectra_raw_timed[i][j]->SetLineColor(colors[i]);
            if (i == 0) {
                spectra_raw_timed[i][j]->Draw();
            } else {
                spectra_raw_timed[i][j]->DrawCopy("same");
            }
            spectra_raw_timed[i][j]->Fit((std::string("gauss_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminbulk[j], fitmaxbulk[j]);
            spectra_raw_timed[i][j]->Fit((std::string("pedest_gauss_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminpedest[j], fitmaxpedest[j]);
            fit_bulk[i][j] = spec_raw_fit[i][j]->GetParameter(1);
            fit_bulk_err[i][j] = spec_raw_fit[i][j]->GetParError(1);
            fit_pedest[i][j] = spec_raw_fit_pedest[i][j]->GetParameter(1);
            fit_pedest_err[i][j] = spec_raw_fit_pedest[i][j]->GetParError(1);

            if (show_fit_bulk || show_fit_peak) {
                spec_raw_fit[i][j]->DrawCopy("same");
                spec_raw_fit_pedest[i][j]->DrawCopy("same");
            }
            if (show_fit_peak) {
                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(fitminpedest[j] - 50, fitmaxpedest[j] + 50);
            } else if (show_fit_bulk) {
                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(fitminbulk[j] - 50, fitmaxbulk[j] + 50);
//                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(400, 800);
                spectra_raw_timed[i][j]->GetYaxis()->SetRangeUser(50, 1000);
            } else {
                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
            }

            gPad->SetLogy();
        }
    }

    can_n_raw_timed->Update();

    // ============================================================================
    // Raw corrected spectra
    // ============================================================================
    auto can_n_corr_timed = new TCanvas("Spectra_corr_n_timed", "Raw corr spectra");

    // Legend area
    TPad *pad_legend_corr = new TPad("pad_legend_corr", "Legend Pad", 0, 0.92, 1, 1);
    pad_legend_corr->SetFillColor(0);
    pad_legend_corr->Draw();
    pad_legend_corr->cd();
    pad_legend_corr->Divide(3, 1);
    std::vector<TLegend*> legend_corr(3);
    pad_legend_corr->cd(1);
    legend_corr[0] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_corr[0]->AddEntry(spectra_corr_timed[0][0], "n channels >= 10#circC", "l");
    legend_corr[0]->AddEntry(spectra_corr_timed[1][0], "n channels 10 - 5#circC", "l");
    legend_corr[0]->Draw();
    pad_legend_corr->cd(2);
    legend_corr[1] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_corr[1]->AddEntry(spectra_corr_timed[2][0], "n channels 5 - 0#circC", "l");
    legend_corr[1]->AddEntry(spectra_corr_timed[3][0], "n channels 0 - -5#circC", "l");
    legend_corr[1]->Draw();
    pad_legend_corr->cd(3);
    legend_corr[2] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_corr[2]->AddEntry(spectra_corr_timed[4][0], "n channels -5 - -10#circC", "l");
    legend_corr[2]->AddEntry(spectra_corr_timed[5][0], "n channels < -10#circC", "l");
    legend_corr[2]->Draw();

    // Histograms area
    can_n_corr_timed->cd();
    TPad *pad_hist_corr = new TPad("pad_hist_corr", "Histogram pad", 0, 0, 1, 0.92);
    pad_hist_corr->Draw();
    pad_hist_corr->cd();
    pad_hist_corr->Divide(8, 4, 0.001, 0.001);

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            pad_hist_corr->cd(j + 1);
            spectra_corr_timed[i][j]->SetStats(0);
            spectra_corr_timed[i][j]->SetLineColor(colors[i]);
            if (i == 0) {
                spectra_corr_timed[i][j]->Draw();
            } else {
                spectra_corr_timed[i][j]->DrawCopy("same");
            }
            spectra_corr_timed[i][j]->Fit((std::string("gauss_corr_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminbulkcorr[j], fitmaxbulkcorr[j]);
            spectra_corr_timed[i][j]->Fit((std::string("pedest_gauss_corr_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminpedestcorr[j], fitmaxpedestcorr[j]);
            fit_corr_bulk[i][j] = spec_corr_fit[i][j]->GetParameter(1);
            fit_corr_bulk_err[i][j] = spec_corr_fit[i][j]->GetParError(1);
            fit_corr_pedest[i][j] = spec_corr_fit_pedest[i][j]->GetParameter(1);
            fit_corr_pedest_err[i][j] = spec_corr_fit_pedest[i][j]->GetParError(1);

            if (show_fit_bulk || show_fit_peak) {
                spec_corr_fit[i][j]->DrawCopy("same");
                spec_corr_fit_pedest[i][j]->DrawCopy("same");
            }
            if (show_fit_peak) {
                spectra_corr_timed[i][j]->GetXaxis()->SetRangeUser(fitminpedestcorr[j] - 50, fitmaxpedestcorr[j] + 50);
            } else if (show_fit_bulk) {
                spectra_corr_timed[i][j]->GetXaxis()->SetRangeUser(fitminbulkcorr[j] - 50, fitmaxbulkcorr[j] + 50);
                spectra_corr_timed[i][j]->GetYaxis()->SetRangeUser(100, 1000);
            } else {
                spectra_corr_timed[i][j]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
            }

            gPad->SetLogy();
        }
    }

    can_n_corr_timed->Update();


    std::vector<std::string> temp_str = {">10", "10-5", "5-0", "0--5", "-5--10", "<-10"};
    for (size_t j=0; j<32; j++) {
        for (size_t i=0; i<6; i++) {
            std::cout << "Bulk mean for temperature : " << temp_str[i] << " channel : " << j << "  -  " << fit_bulk[i][j] << "+-" << fit_bulk_err[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    for (size_t j=0; j<32; j++) {
        for (size_t i=0; i<6; i++) {
            std::cout << "Peak mean for temperature : " << temp_str[i] << " channel : " << j << "  -  " << fit_pedest[i][j] << "+-" << fit_pedest_err[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    // Open the bulk and peak fit files
    std::ofstream bulk_file("./Kiruna_data/calib_dssd/fit_bulk.txt");
    if (bulk_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                bulk_file << fit_bulk[i][j];
                if (i < 5) {
                    bulk_file << " ";
                } else {
                    bulk_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                bulk_file << fit_bulk_err[i][j];
                if (i < 5) {
                    bulk_file << " ";
                }
            }

            bulk_file << "\n";
        }
        bulk_file.close();
        std::cout << "Bulk fit data saved in fit_bulk.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }

    std::ofstream peak_file("./Kiruna_data/calib_dssd/fit_pedest.txt");
    if (peak_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                peak_file << fit_pedest[i][j];
                if (i < 5) {
                    peak_file << " ";
                } else {
                    peak_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                peak_file << fit_pedest_err[i][j];
                if (i < 5) {
                    peak_file << " ";
                }
            }

            peak_file << "\n";
        }
        peak_file.close();
        std::cout << "Peak fit data saved in fit_pedest.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }

    // Open the corr_bulk and corr_peak fit files
    std::ofstream bulk_corr_file("./Kiruna_data/calib_dssd/fit_corr_bulk.txt");
    if (bulk_corr_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                bulk_corr_file << fit_corr_bulk[i][j];
                if (i < 5) {
                    bulk_corr_file << " ";
                } else {
                    bulk_corr_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                bulk_corr_file << fit_corr_bulk_err[i][j];
                if (i < 5) {
                    bulk_corr_file << " ";
                }
            }

            bulk_corr_file << "\n";
        }
        bulk_corr_file.close();
        std::cout << "Bulk corr fit data saved in fit_corr_bulk.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }

    std::ofstream peak_corr_file("./Kiruna_data/calib_dssd/fit_corr_pedest.txt");
    if (peak_corr_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                peak_corr_file << fit_corr_pedest[i][j];
                if (i < 5) {
                    peak_corr_file << " ";
                } else {
                    peak_corr_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                peak_corr_file << fit_corr_pedest_err[i][j];
                if (i < 5) {
                    peak_corr_file << " ";
                }
            }

            peak_corr_file << "\n";
        }
        peak_corr_file.close();
        std::cout << "Peak corr fit data saved in fit_corr_pedest.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }


    // ============================================================================
    // For the manuscript !
    // ============================================================================
    // ============================================================================
    // Calibrated spectra
    // ============================================================================
    auto can_chan13_calib = new TCanvas("Spectra_chan13_calibrated", "Spectra 13 calibrated", 900, 600);

    // Legend area
    TPad *pad_legend_calib_13 = new TPad("pad_legend_calib", "Legend Pad", 0, 0.9, 1, 1);
    pad_legend_calib_13->SetFillColor(0);
    pad_legend_calib_13->Draw();
    pad_legend_calib_13->cd();
    pad_legend_calib_13->Divide(3, 1);
    pad_legend_calib_13->cd(1);
    legend_calib[0]->SetTextSize(0.33);
    legend_calib[0]->Draw();
    pad_legend_calib_13->cd(2);
    legend_calib[1]->SetTextSize(0.33);
    legend_calib[1]->Draw();
    pad_legend_calib_13->cd(3);
    legend_calib[2]->SetTextSize(0.33);
    legend_calib[2]->Draw();

    // Histograms area
    can_chan13_calib->cd();
    TPad *pad_hist13_calib = new TPad("pad_hist13_calib", "Histogram 13 pad", 0, 0, 1, 0.9);
    pad_hist13_calib->Draw();
    pad_hist13_calib->cd();

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        spectra_calib_timed[i][13]->SetStats(0);
        spectra_calib_timed[i][13]->GetXaxis()->SetRangeUser(xlimmincalib, xlimmaxcalib);
        spectra_calib_timed[i][13]->SetLineColor(colors[i]);
        spectra_calib_timed[i][13]->GetXaxis()->SetTitleSize(0.045);
        spectra_calib_timed[i][13]->GetYaxis()->SetTitleSize(0.045);
        spectra_calib_timed[i][13]->GetXaxis()->SetLabelSize(0.04);
        spectra_calib_timed[i][13]->GetYaxis()->SetLabelSize(0.04);
        spectra_calib_timed[i][13]->GetXaxis()->SetTitle("Energy (keV)");
        spectra_calib_timed[i][13]->GetYaxis()->SetTitle("Number of events");
        spectra_calib_timed[i][13]->SetTitle("");
        gPad->SetTopMargin(0.05);
//            spectra_calib_timed[i][13]->Fit((std::string("gauss_fit_timed1") + std::to_string(i)).c_str(), "", "", fitxmin2, fitxmax2);
        if (i == 0) {
            spectra_calib_timed[i][13]->Draw();
        } else {
            spectra_calib_timed[i][13]->DrawCopy("same");
        }
        gPad->SetLogy();
    }

    can_chan13_calib->Update();

    // ============================================================================
    // Raw spectra
    // ============================================================================
    auto can_chan13_raw = new TCanvas("Spectra_chan13_raw", "Raw spectra chan13", 900, 600);

    // Legend area
    TPad *pad_legend_raw_13 = new TPad("pad_legend_raw13", "Legend Pad13", 0, 0.9, 1, 1);
    pad_legend_raw_13->SetFillColor(0);
    pad_legend_raw_13->Draw();
    pad_legend_raw_13->cd();
    pad_legend_raw_13->Divide(3, 1);
    pad_legend_raw_13->cd(1);
    legend_raw[0]->SetTextSize(0.33);
    legend_raw[0]->Draw();
    pad_legend_raw_13->cd(2);
    legend_raw[1]->SetTextSize(0.33);
    legend_raw[1]->Draw();
    pad_legend_raw_13->cd(3);
    legend_raw[2]->SetTextSize(0.33);
    legend_raw[2]->Draw();

    // Histograms area
    can_chan13_raw->cd();
    TPad *pad_hist13_raw = new TPad("pad_hist13_raw", "Histogram 13 pad", 0, 0, 1, 0.9);
    pad_hist13_raw->Draw();
    pad_hist13_raw->cd();

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        spectra_raw_timed[i][13]->SetStats(0);
        spectra_raw_timed[i][13]->SetLineColor(colors[i]);
        spectra_raw_timed[i][13]->GetXaxis()->SetTitleSize(0.045);
        spectra_raw_timed[i][13]->GetYaxis()->SetTitleSize(0.045);
        spectra_raw_timed[i][13]->GetXaxis()->SetLabelSize(0.04);
        spectra_raw_timed[i][13]->GetYaxis()->SetLabelSize(0.04);
        spectra_raw_timed[i][13]->GetXaxis()->SetTitle("Energy (ADC)");
        spectra_raw_timed[i][13]->GetYaxis()->SetTitle("Number of events");
        spectra_raw_timed[i][13]->SetTitle("");
        gPad->SetTopMargin(0.05);
        if (i == 0) {
            spectra_raw_timed[i][13]->Draw();
        } else {
            spectra_raw_timed[i][13]->DrawCopy("same");
        }
        spectra_raw_timed[i][13]->Fit((std::string("gauss_fit") + std::to_string(i) + std::string("_13")).c_str(), "QNR+", "", fitminbulk[13], fitmaxbulk[13]);
        spectra_raw_timed[i][13]->Fit((std::string("pedest_gauss_fit") + std::to_string(i) + std::string("_13")).c_str(), "QNR+", "", fitminpedest[13], fitmaxpedest[13]);
        fit_bulk[i][13] = spec_raw_fit[i][13]->GetParameter(1);
        fit_bulk_err[i][13] = spec_raw_fit[i][13]->GetParError(1);
        fit_pedest[i][13] = spec_raw_fit_pedest[i][13]->GetParameter(1);
        fit_pedest_err[i][13] = spec_raw_fit_pedest[i][13]->GetParError(1);

//        if (show_fit_bulk || show_fit_peak) {
    //        spec_raw_fit[i][13]->DrawCopy("same");
    //        spec_raw_fit_pedest[i][13]->DrawCopy("same");
//        }
//        if (show_fit_peak) {
//            spectra_raw_timed[i][13]->GetXaxis()->SetRangeUser(fitminpedest[13] - 50, fitmaxpedest[13] + 50);
//        } else if (show_fit_bulk) {
//            spectra_raw_timed[i][13]->GetXaxis()->SetRangeUser(fitminbulk[13] - 50, fitmaxbulk[13] + 50);
//                spectra_raw_timed[i][13]->GetXaxis()->SetRangeUser(400, 800);
//            spectra_raw_timed[i][13]->GetYaxis()->SetRangeUser(50, 1000);
//        } else {
//            spectra_raw_timed[i][13]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
//        }

        gPad->SetLogy();

    }

    can_chan13_raw->Update();





////  =============== 2D spectra ===============
//    std::cout << "Filling 2D spectra" << std::endl;
//    auto hist2D1 = new TH2D("hist2D1", "ADC channel for time vs strip number map", 64, 0, 64, end_time - init_time, init_time, end_time);
//    auto hist2D2 = new TH2D("hist2D2", "adc channel vs strip number", 64, 0, 64, 1024, 0, 1024);
//    auto hist2D1calib = new TH2D("hist2D1calib", "energy channel for time vs strip number map", 64, 0, 64, end_time - init_time, init_time, end_time);
//    auto hist2D2calib = new TH2D("hist2D2calib", "energy channel vs strip number", 64, 0, 64, 1500, 0, 1500);
//
//    for (uint64_t i = 0; i < nentries; ++i) {
//        tree_hist->GetEntry(i);
////        std::cout << "Entry : " << time_corr_abs << " OK " << std::endl;
//        for (uint64_t ch = 0; ch < 64; ++ch) {
////            std::cout << adc_channels[ch] << "  -  " << adc_channels_calib[ch] << std::endl;
//            hist2D1->Fill(static_cast<Double_t>(ch), static_cast<Double_t>(i), static_cast<Double_t>(adc_channels[ch]));
//            hist2D2->Fill(ch, adc_channels[ch]);
//            hist2D1calib->Fill(static_cast<Double_t>(ch), static_cast<Double_t>(i), static_cast<Double_t>(adc_channels_calib[ch]));
//            if (ch != 32 && ch != 34 && ch != 38 && ch != 45 && ch != 56 && ch != 63) {
//                hist2D2calib->Fill(ch, adc_channels_calib[ch]);
//            } else {
//                hist2D2calib->Fill(ch, 0);
//            }
//        }
//    }
//
//    auto can1 = new TCanvas("canvas", "DSSD channels 2D hist");
//    can1->Divide(2, 2);
//    can1->cd(1);
//    hist2D1->SetStats(0);
//    hist2D1->Draw("COLZ");
//
//    can1->cd(2);
//    hist2D2->GetXaxis()->SetRangeUser(32, 64);
//    hist2D2->GetYaxis()->SetRangeUser(420, 1024);
//    hist2D2->SetStats(0);
//    hist2D2->Draw("COLZ");
//
//    can1->cd(3);
//    hist2D1calib->SetStats(0);
//    hist2D1calib->Draw("COLZ");
//
//    can1->cd(4);
//    gPad->SetLogz();
//    hist2D2calib->GetXaxis()->SetRangeUser(32, 64);
//    hist2D2calib->GetYaxis()->SetRangeUser(xlimmin, 1200);
//    hist2D2calib->SetStats(0);
//    hist2D2calib->Draw("COLZ");

//    can1->SaveAs("DSSD_channels_2D_hist.png");

    DisplayApp.Run();

}

void DisplayDSSDChannels2(const std::string &tree_name) {
    // Extraction of the root tree
    file_exists(tree_name);
    auto tree_file = new TFile(tree_name.c_str());
    if (!tree_file || tree_file->IsZombie()) {
        std::cerr << "Error while opening the file." << std::endl;
    }
    auto tree_hist = (TTree*) tree_file->Get("Events");

    // Variables to read the initial tree
    Double_t time_corr_abs;
    uint16_t adc_channels[64];
    Double_t adc_channels_calib[64];
    Double_t adc_channels_calib_corr[64];
    Double_t temperature;
    tree_hist->SetBranchAddress("time_corr_abs",     &time_corr_abs);
    tree_hist->SetBranchAddress("adc_channels",      &adc_channels);
    tree_hist->SetBranchAddress("adc_channels_calib",      &adc_channels_calib);
    tree_hist->SetBranchAddress("adc_channels_calib_corr",      &adc_channels_calib_corr);
    tree_hist->SetBranchAddress("temperature", &temperature);

    uint64_t nentries = tree_hist->GetEntries();
    uint32_t init_time;
    uint32_t end_time;
    tree_hist->GetEntry(0);
    init_time = static_cast<uint32_t>(time_corr_abs);
    tree_hist->GetEntry(nentries - 1);
    end_time = static_cast<uint32_t>(time_corr_abs + 1) ;


    TApplication DisplayApp("App", nullptr, nullptr);

    int xlimmin = 0;
    int xlimmax = 1500;
    int xlimmincalib = 0;
    int xlimmaxcalib = 1500;

//  =============== time vs temperature ===============
//    auto can_time_temp = new TCanvas("time_vs_temp", "Time vs temperature");
//    tree_hist->Draw("temperature:time_corr_abs");

    bool show_fit_bump = false;
    Double_t fitxmin2 = 450;
    Double_t fitxmax2 = 700;

    Double_t fit_bump_simple_min = 400;
    Double_t fit_bump_simple_max = 800;
    std::vector<Double_t> fitminbumpraw = {420, 460, 400, 400, 450, 430, 420, 420, 450, 430, 410, 420, 420, 440, 430, 440, 440, 420, 420, 420, 420, 420, 430, 430, 420, 410, 440, 430, 430, 420, 430, 440};
    std::vector<Double_t> fitmaxbumpraw = {750, 800, 800, 770, 790, 780, 750, 760, 770, 770, 760, 780, 760, 770, 750, 750, 760, 770, 770, 770, 770, 760, 760, 760, 750, 740, 770, 760, 760, 750, 760, 750};
    std::vector<Double_t> fitminbumpcalib = {420, 460, 520, 430, 440, 420, 420, 430, 440, 430, 420, 420, 420, 420, 430, 430, 430, 420, 420, 420, 420, 420, 420, 420, 430, 430, 430, 420, 430, 430, 440, 440};
    std::vector<Double_t> fitmaxbumpcalib = {720, 760, 820, 730, 740, 720, 720, 730, 740, 730, 720, 720, 720, 720, 730, 730, 730, 720, 720, 720, 720, 720, 720, 720, 730, 730, 730, 720, 730, 730, 740, 740};
    std::vector<Double_t> fitminbumpcalibcorr = {480, 500, 490, 480, 490, 480, 480, 470,  480, 480, 480, 480, 480, 480, 470, 480,  480, 470, 470, 480, 480, 480, 480, 480,  470, 470, 480, 480, 480, 480, 500, 480};
    std::vector<Double_t> fitmaxbumpcalibcorr = {700, 720, 710, 690, 700, 690, 690, 690,  690, 690, 690, 690, 690, 690, 690, 690,  690, 690, 690, 690, 690, 690, 690, 690,  690, 690, 700, 690, 690, 690, 710, 690};

    // Spec histograms
    std::vector<std::vector<TH1D*>> spectra_raw_timed(6, std::vector<TH1D*>(32));
    std::vector<std::vector<TH1D*>> spectra_calib_timed(6, std::vector<TH1D*>(32));
    std::vector<std::vector<TH1D*>> spectra_calib_corr_timed(6, std::vector<TH1D*>(32));
    // Spec fits
    std::vector<std::vector<TF1*>> spec_raw_fit(6, std::vector<TF1*>(32));
    std::vector<std::vector<TF1*>> spec_calib_fit(6, std::vector<TF1*>(32));
    std::vector<std::vector<TF1*>> spec_calib_corr_fit(6, std::vector<TF1*>(32));
    // Fit calib channels
    std::vector<std::vector<Double_t>> fit_raw_bump(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_raw_bump_err(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_calib_bump(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_calib_bump_err(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_calib_corr_bump(6, std::vector<Double_t>(32));
    std::vector<std::vector<Double_t>> fit_calib_corr_bump_err(6, std::vector<Double_t>(32));
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            // Spec histograms
            spectra_raw_timed[i][j] = new TH1D((std::string("spec_raw_timed") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "Raw ADC spectrum", 1024, 0, 1024);
            spectra_calib_timed[i][j] = new TH1D((std::string("spec_calib_timed") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "Energy spectrum", 250, 0, 1500);
            spectra_calib_corr_timed[i][j] = new TH1D((std::string("spec_calib_corr_timed") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "Corr energy spectrum", 250, 0, 1500);
            // Initializing the calib fits
            spec_raw_fit[i][j] = new TF1((std::string("gauss_raw_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminbumpraw[j], fitmaxbumpraw[j]);
            spec_raw_fit[i][j]->SetLineColor(kRed);
            // Initializing the calib fits
            spec_calib_fit[i][j] = new TF1((std::string("gauss_calib_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminbumpcalib[j], fitmaxbumpcalib[j]);
            spec_calib_fit[i][j]->SetLineColor(kRed);
            // Initializing the calib corr fits
            spec_calib_corr_fit[i][j] = new TF1((std::string("gauss_calib_corr_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "gaus", fitminbumpcalibcorr[j], fitmaxbumpcalibcorr[j]);
            spec_calib_corr_fit[i][j]->SetLineColor(kRed);
        }
    }

    // Filling histograms
    size_t temp_ite;
    for (uint64_t i = 0; i < nentries; ++i) {
        tree_hist->GetEntry(i);
        if (temperature >= 10) {
            temp_ite = 0;
        } else if (temperature < 10 && temperature >= 5) {
            temp_ite = 1;
        } else if (temperature < 5 && temperature >= 0) {
            temp_ite = 2;
        } else if (temperature < 0 && temperature >= -5) {
            temp_ite = 3;
        } else if (temperature < -5 && temperature >= -10) {
            temp_ite = 4;
        } else if (temperature < -10) {
            temp_ite = 5;
        }
        for (size_t j=0; j<32; j++) {
            spectra_raw_timed[temp_ite][j]->Fill(adc_channels[j + 32]);
            spectra_calib_timed[temp_ite][j]->Fill(adc_channels_calib[j + 32]);
            spectra_calib_corr_timed[temp_ite][j]->Fill(adc_channels_calib_corr[j + 32]);
        }
    }

    // Colors for different energy ranges
    int colors[6] = {kRed, kBlue, kGreen+2, kMagenta, kOrange, kCyan};

    // ============================================================================
    // Raw spectra
    // ============================================================================
    auto can_n_raw_timed = new TCanvas("Spectra_raw_n_timed", "Raw spectra");

    // Legend area
    TPad *pad_legend_raw = new TPad("pad_legend_raw", "Legend Pad", 0, 0.92, 1, 1);
    pad_legend_raw->SetFillColor(0);
    pad_legend_raw->Draw();
    pad_legend_raw->cd();
    pad_legend_raw->Divide(3, 1);
    std::vector<TLegend*> legend_raw(3);
    pad_legend_raw->cd(1);
    legend_raw[0] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_raw[0]->AddEntry(spectra_raw_timed[0][0], "n channels >= 10#circC", "l");
    legend_raw[0]->AddEntry(spectra_raw_timed[1][0], "n channels 10 - 5#circC", "l");
    legend_raw[0]->Draw();
    pad_legend_raw->cd(2);
    legend_raw[1] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_raw[1]->AddEntry(spectra_raw_timed[2][0], "n channels 5 - 0#circC", "l");
    legend_raw[1]->AddEntry(spectra_raw_timed[3][0], "n channels 0 - -5#circC", "l");
    legend_raw[1]->Draw();
    pad_legend_raw->cd(3);
    legend_raw[2] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_raw[2]->AddEntry(spectra_raw_timed[4][0], "n channels -5 - -10#circC", "l");
    legend_raw[2]->AddEntry(spectra_raw_timed[5][0], "n channels < -10#circC", "l");
    legend_raw[2]->Draw();

    // Histograms area
    can_n_raw_timed->cd();
    TPad *pad_hist_raw = new TPad("pad_hist_raw", "Histogram pad", 0, 0, 1, 0.92);
    pad_hist_raw->Draw();
    pad_hist_raw->cd();
    pad_hist_raw->Divide(8, 4, 0.001, 0.001);

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            pad_hist_raw->cd(j + 1);
            spectra_raw_timed[i][j]->SetStats(0);
            spectra_raw_timed[i][j]->SetLineColor(colors[i]);
            if (i == 0) {
                spectra_raw_timed[i][j]->Draw();
            } else {
                spectra_raw_timed[i][j]->DrawCopy("same");
            }

            spectra_raw_timed[i][j]->Fit((std::string("gauss_raw_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminbumpraw[j], fitmaxbumpraw[j]);
            fit_raw_bump[i][j] = spec_raw_fit[i][j]->GetParameter(1);
            fit_raw_bump_err[i][j] = spec_raw_fit[i][j]->GetParError(1);
            if (show_fit_bump) {
                spec_raw_fit[i][j]->DrawCopy("same");
                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(fitminbumpraw[j] - 50, fitmaxbumpraw[j] + 50);
                spectra_raw_timed[i][j]->GetYaxis()->SetRangeUser(50, 1000);
            } else {
//                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
                spectra_raw_timed[i][j]->GetXaxis()->SetRangeUser(fitminbumpraw[j], fitmaxbumpraw[j]);
            }

            gPad->SetLogy();
        }
    }
    can_n_raw_timed->Update();

    // ============================================================================
    // Calibrated spectra
    // ============================================================================
    auto can_n_calib_timed = new TCanvas("Spectra_calibrated_n_timed", "Spectra calibrated");

    // Legend area
    TPad *pad_legend_calib = new TPad("pad_legend_calib", "Legend Pad", 0, 0.92, 1, 1);
    pad_legend_calib->SetFillColor(0);
    pad_legend_calib->Draw();
    pad_legend_calib->cd();
    pad_legend_calib->Divide(3, 1);
    std::vector<TLegend*> legend_calib(3);
    pad_legend_calib->cd(1);
    legend_calib[0] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_calib[0]->AddEntry(spectra_calib_timed[0][0], "n channels >= 10#circC", "l");
    legend_calib[0]->AddEntry(spectra_calib_timed[1][0], "n channels 10 - 5#circC", "l");
    legend_calib[0]->Draw();
    pad_legend_calib->cd(2);
    legend_calib[1] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_calib[1]->AddEntry(spectra_calib_timed[2][0], "n channels 5 - 0#circC", "l");
    legend_calib[1]->AddEntry(spectra_calib_timed[3][0], "n channels 0 - -5#circC", "l");
    legend_calib[1]->Draw();
    pad_legend_calib->cd(3);
    legend_calib[2] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_calib[2]->AddEntry(spectra_calib_timed[4][0], "n channels -5 - -10#circC", "l");
    legend_calib[2]->AddEntry(spectra_calib_timed[5][0], "n channels < -10#circC", "l");
    legend_calib[2]->Draw();

    // Histograms area
    can_n_calib_timed->cd();
    TPad *pad_hist_calib = new TPad("pad_hist_calib", "Histogram pad", 0, 0, 1, 0.92);
    pad_hist_calib->Draw();
    pad_hist_calib->cd();
    pad_hist_calib->Divide(8, 4, 0.001, 0.001);

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            pad_hist_calib->cd(j + 1);
            spectra_calib_timed[i][j]->SetStats(0);
            spectra_calib_timed[i][j]->GetXaxis()->SetRangeUser(xlimmincalib, xlimmaxcalib);
            spectra_calib_timed[i][j]->SetLineColor(colors[i]);
            if (i == 0) {
                spectra_calib_timed[i][j]->Draw();
            } else {
                spectra_calib_timed[i][j]->DrawCopy("same");
            }

            spectra_calib_timed[i][j]->Fit((std::string("gauss_calib_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminbumpcalib[j], fitmaxbumpcalib[j]);
            fit_calib_bump[i][j] = spec_calib_fit[i][j]->GetParameter(1);
            fit_calib_bump_err[i][j] = spec_calib_fit[i][j]->GetParError(1);
            if (show_fit_bump) {
                spec_calib_fit[i][j]->DrawCopy("same");
                spectra_calib_timed[i][j]->GetXaxis()->SetRangeUser(fitminbumpraw[j] - 50, fitmaxbumpraw[j] + 50);
                spectra_calib_timed[i][j]->GetYaxis()->SetRangeUser(50, 1000);
            } else {
//                spectra_calib_timed[i][j]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
                spectra_calib_timed[i][j]->GetXaxis()->SetRangeUser(fitminbumpcalib[j], fitmaxbumpcalib[j]);
            }

            gPad->SetLogy();
        }
    }

    can_n_calib_timed->Update();

    // ============================================================================
    // Calib corrected spectra
    // ============================================================================
    auto can_n_corr_timed = new TCanvas("Spectra_corr_n_timed", "Calib corr spectra");

    // Legend area
    TPad *pad_legend_corr = new TPad("pad_legend_corr", "Legend Pad", 0, 0.92, 1, 1);
    pad_legend_corr->SetFillColor(0);
    pad_legend_corr->Draw();
    pad_legend_corr->cd();
    pad_legend_corr->Divide(3, 1);
    std::vector<TLegend*> legend_corr(3);
    pad_legend_corr->cd(1);
    legend_corr[0] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_corr[0]->AddEntry(spectra_calib_corr_timed[0][0], "n channels >= 10#circC", "l");
    legend_corr[0]->AddEntry(spectra_calib_corr_timed[1][0], "n channels 10 - 5#circC", "l");
    legend_corr[0]->Draw();
    pad_legend_corr->cd(2);
    legend_corr[1] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_corr[1]->AddEntry(spectra_calib_corr_timed[2][0], "n channels 5 - 0#circC", "l");
    legend_corr[1]->AddEntry(spectra_calib_corr_timed[3][0], "n channels 0 - -5#circC", "l");
    legend_corr[1]->Draw();
    pad_legend_corr->cd(3);
    legend_corr[2] = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend_corr[2]->AddEntry(spectra_calib_corr_timed[4][0], "n channels -5 - -10#circC", "l");
    legend_corr[2]->AddEntry(spectra_calib_corr_timed[5][0], "n channels < -10#circC", "l");
    legend_corr[2]->Draw();

    // Histograms area
    can_n_corr_timed->cd();
    TPad *pad_hist_corr = new TPad("pad_hist_corr", "Histogram pad", 0, 0, 1, 0.92);
    pad_hist_corr->Draw();
    pad_hist_corr->cd();
    pad_hist_corr->Divide(8, 4, 0.001, 0.001);

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<32; j++) {
            pad_hist_corr->cd(j + 1);
            spectra_calib_corr_timed[i][j]->SetStats(0);
            spectra_calib_corr_timed[i][j]->SetLineColor(colors[i]);
            if (i == 0) {
                spectra_calib_corr_timed[i][j]->Draw();
            } else {
                spectra_calib_corr_timed[i][j]->DrawCopy("same");
            }
            spectra_calib_corr_timed[i][j]->Fit((std::string("gauss_calib_corr_fit") + std::to_string(i) + std::string("_") + std::to_string(j)).c_str(), "QNR+", "", fitminbumpcalibcorr[j], fitmaxbumpcalibcorr[j]);
            fit_calib_corr_bump[i][j] = spec_calib_corr_fit[i][j]->GetParameter(1);
            fit_calib_corr_bump_err[i][j] = spec_calib_corr_fit[i][j]->GetParError(1);

            if (show_fit_bump) {
                spec_calib_corr_fit[i][j]->DrawCopy("same");
                spectra_calib_corr_timed[i][j]->GetXaxis()->SetRangeUser(fitminbumpcalibcorr[j] - 50, fitminbumpcalibcorr[j] + 50);
                spectra_calib_corr_timed[i][j]->GetYaxis()->SetRangeUser(100, 1000);
            } else {
                spectra_calib_corr_timed[i][j]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
            }

            gPad->SetLogy();
        }
    }
    can_n_corr_timed->Update();


    std::vector<std::string> temp_str = {">10", "10-5", "5-0", "0--5", "-5--10", "<-10"};
    for (size_t j=0; j<32; j++) {
        for (size_t i=0; i<6; i++) {
            std::cout << "Bump raw mean, channel : " << j << "  -  " << fit_raw_bump[i][j] << "+-" << fit_raw_bump_err[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    // Open the bump and peak fit files
    std::ofstream bump_file("./Kiruna_data/calib_dssd/fit_raw_bump.txt");
    if (bump_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                bump_file << fit_raw_bump[i][j];
                if (i < 5) {
                    bump_file << " ";
                } else {
                    bump_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                bump_file << fit_raw_bump_err[i][j];
                if (i < 5) {
                    bump_file << " ";
                }
            }

            bump_file << "\n";
        }
        bump_file.close();
        std::cout << "Bump fit data saved in fit_bump.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }

    for (size_t j=0; j<32; j++) {
        for (size_t i=0; i<6; i++) {
            std::cout << "Bump calib at 25°C, channel : " << j << "  -  " << fit_calib_bump[i][j] << "+-" << fit_calib_bump_err[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    // Open the bump and peak fit files
    std::ofstream bump_calib_file("./Kiruna_data/calib_dssd/fit_calib_bump.txt");
    if (bump_calib_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                bump_calib_file << fit_calib_bump[i][j];
                if (i < 5) {
                    bump_calib_file << " ";
                } else {
                    bump_calib_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                bump_calib_file << fit_calib_bump_err[i][j];
                if (i < 5) {
                    bump_calib_file << " ";
                }
            }

            bump_calib_file << "\n";
        }
        bump_calib_file.close();
        std::cout << "Bump fit data saved in fit_calib_bump.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }

    for (size_t j=0; j<32; j++) {
        for (size_t i=0; i<6; i++) {
            std::cout << "Bump mean at detection temperature : " << temp_str[i] << " channel : " << j << "  -  " << fit_calib_corr_bump[i][j] << "+-" << fit_calib_corr_bump_err[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    // Open the corr_bump and corr_peak fit files
    std::ofstream bump_corr_file("./Kiruna_data/calib_dssd/fit_calib_corr_bump.txt");
    if (bump_corr_file.is_open()) {
        for (size_t j = 0; j < 32; j++) {
            for (size_t i = 0; i < 6; i++) {
                bump_corr_file << fit_calib_corr_bump[i][j];
                if (i < 5) {
                    bump_corr_file << " ";
                } else {
                    bump_corr_file << "|";
                }
            }
            for (size_t i = 0; i < 6; i++) {
                bump_corr_file << fit_calib_corr_bump_err[i][j];
                if (i < 5) {
                    bump_corr_file << " ";
                }
            }

            bump_corr_file << "\n";
        }
        bump_corr_file.close();
        std::cout << "Bump corr fit data saved in fit_corr_bump.txt\n";
    } else {
        std::cerr << "Error : file not opening !\n";
    }


    // ============================================================================
    // For the manuscript !
    // ============================================================================
    size_t specite = 20;

    // ============================================================================
    // Raw spectra
    // ============================================================================
    auto can_chani_raw = new TCanvas((std::string("Spectra_chan") + std::to_string(specite) + std::string("_raw")).c_str(), "Raw spectrum", 900, 600);

    // Legend area
    TPad *pad_legend_raw_i = new TPad("pad_legend_rawi", "Legend Padi", 0, 0.9, 1, 1);
    pad_legend_raw_i->SetFillColor(0);
    pad_legend_raw_i->Draw();
    pad_legend_raw_i->cd();
    pad_legend_raw_i->Divide(3, 1);
    pad_legend_raw_i->cd(1);
    legend_raw[0]->SetTextSize(0.33);
    legend_raw[0]->Draw();
    pad_legend_raw_i->cd(2);
    legend_raw[1]->SetTextSize(0.33);
    legend_raw[1]->Draw();
    pad_legend_raw_i->cd(3);
    legend_raw[2]->SetTextSize(0.33);
    legend_raw[2]->Draw();

    // Histograms area
    can_chani_raw->cd();
    TPad *pad_histi_raw = new TPad("pad_histi_raw", "Histogram pad", 0, 0, 1, 0.9);
    pad_histi_raw->Draw();
    pad_histi_raw->cd();

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        spectra_raw_timed[i][specite]->SetStats(0);
        spectra_raw_timed[i][specite]->GetXaxis()->SetRangeUser(0, 1024);
        spectra_raw_timed[i][specite]->SetLineColor(colors[i]);
        spectra_raw_timed[i][specite]->GetXaxis()->SetTitleSize(0.045);
        spectra_raw_timed[i][specite]->GetYaxis()->SetTitleSize(0.045);
        spectra_raw_timed[i][specite]->GetXaxis()->SetLabelSize(0.04);
        spectra_raw_timed[i][specite]->GetYaxis()->SetLabelSize(0.04);
        spectra_raw_timed[i][specite]->GetXaxis()->SetTitle("Energy (ADC)");
        spectra_raw_timed[i][specite]->GetYaxis()->SetTitle("Number of events");
        spectra_raw_timed[i][specite]->SetTitle("");
        gPad->SetTopMargin(0.05);

        if (i == 0) {
            spectra_raw_timed[i][specite]->Draw();
        } else {
            spectra_raw_timed[i][specite]->DrawCopy("same");
        }
        spectra_raw_timed[i][specite]->Fit((std::string("gauss_fit") + std::to_string(i) + std::string("_") + std::to_string(specite)).c_str(), "QNR+", "", fitminbumpraw[specite], fitmaxbumpraw[specite]);
//        fit_bump[i][specite] = spec_raw_fit[i][specite]->GetParameter(1);
//        fit_bump_err[i][specite] = spec_raw_fit[i][specite]->GetParError(1);
        gPad->SetLogy();
    }
    can_chani_raw->Update();

    // ============================================================================
    // Calibrated spectra
    // ============================================================================
    auto can_chani_calib = new TCanvas((std::string("Spectra_chan") + std::to_string(specite) + std::string("_calib")).c_str(), "Spectrum calibrated", 900, 600);

    // Legend area
    TPad *pad_legend_calib_i = new TPad("pad_legend_calib", "Legend Pad", 0, 0.9, 1, 1);
    pad_legend_calib_i->SetFillColor(0);
    pad_legend_calib_i->Draw();
    pad_legend_calib_i->cd();
    pad_legend_calib_i->Divide(3, 1);
    pad_legend_calib_i->cd(1);
    legend_calib[0]->SetTextSize(0.33);
    legend_calib[0]->Draw();
    pad_legend_calib_i->cd(2);
    legend_calib[1]->SetTextSize(0.33);
    legend_calib[1]->Draw();
    pad_legend_calib_i->cd(3);
    legend_calib[2]->SetTextSize(0.33);
    legend_calib[2]->Draw();

    // Histograms area
    can_chani_calib->cd();
    TPad *pad_histi_calib = new TPad("pad_histi_calib", "Histogram i pad", 0, 0, 1, 0.9);
    pad_histi_calib->Draw();
    pad_histi_calib->cd();

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        spectra_calib_timed[i][specite]->SetStats(0);
        spectra_calib_timed[i][specite]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
        spectra_calib_timed[i][specite]->SetLineColor(colors[i]);
        spectra_calib_timed[i][specite]->GetXaxis()->SetTitleSize(0.045);
        spectra_calib_timed[i][specite]->GetYaxis()->SetTitleSize(0.045);
        spectra_calib_timed[i][specite]->GetXaxis()->SetLabelSize(0.04);
        spectra_calib_timed[i][specite]->GetYaxis()->SetLabelSize(0.04);
        spectra_calib_timed[i][specite]->GetXaxis()->SetTitle("Energy (keV)");
        spectra_calib_timed[i][specite]->GetYaxis()->SetTitle("Number of events");
        spectra_calib_timed[i][specite]->SetTitle("");
        gPad->SetTopMargin(0.05);
        if (i == 0) {
            spectra_calib_timed[i][specite]->Draw();
        } else {
            spectra_calib_timed[i][specite]->DrawCopy("same");
        }
        gPad->SetLogy();
    }
    can_chani_calib->Update();

    // ============================================================================
    // Calibrated corrected spectra
    // ============================================================================
    auto can_chani_calibcorr = new TCanvas((std::string("Spectra_chan") + std::to_string(specite) + std::string("_calibcorr")).c_str(), "Spectrum calibrated corr", 900, 600);

    // Legend area
    TPad *pad_legend_calibcorr_i = new TPad("pad_legend_calibcorr", "Legend Pad", 0, 0.9, 1, 1);
    pad_legend_calibcorr_i->SetFillColor(0);
    pad_legend_calibcorr_i->Draw();
    pad_legend_calibcorr_i->cd();
    pad_legend_calibcorr_i->Divide(3, 1);
    pad_legend_calibcorr_i->cd(1);
    legend_corr[0]->SetTextSize(0.33);
    legend_corr[0]->Draw();
    pad_legend_calibcorr_i->cd(2);
    legend_corr[1]->SetTextSize(0.33);
    legend_corr[1]->Draw();
    pad_legend_calibcorr_i->cd(3);
    legend_corr[2]->SetTextSize(0.33);
    legend_corr[2]->Draw();

    // Histograms area
    can_chani_calibcorr->cd();
    TPad *pad_histi_calibcorr = new TPad("pad_histi_calibcorr", "Histogram i pad", 0, 0, 1, 0.9);
    pad_histi_calibcorr->Draw();
    pad_histi_calibcorr->cd();

    // Drawing the different spectra on different pads
    for (size_t i=0; i<6; i++) {
        spectra_calib_corr_timed[i][specite]->SetStats(0);
        spectra_calib_corr_timed[i][specite]->GetXaxis()->SetRangeUser(xlimmin, xlimmax);
        spectra_calib_corr_timed[i][specite]->SetLineColor(colors[i]);
        spectra_calib_corr_timed[i][specite]->GetXaxis()->SetTitleSize(0.045);
        spectra_calib_corr_timed[i][specite]->GetYaxis()->SetTitleSize(0.045);
        spectra_calib_corr_timed[i][specite]->GetXaxis()->SetLabelSize(0.04);
        spectra_calib_corr_timed[i][specite]->GetYaxis()->SetLabelSize(0.04);
        spectra_calib_corr_timed[i][specite]->GetXaxis()->SetTitle("Energy (keV)");
        spectra_calib_corr_timed[i][specite]->GetYaxis()->SetTitle("Number of events");
        spectra_calib_corr_timed[i][specite]->SetTitle("");
        gPad->SetTopMargin(0.05);
        if (i == 0) {
            spectra_calib_corr_timed[i][specite]->Draw();
        } else {
            spectra_calib_corr_timed[i][specite]->DrawCopy("same");
        }
        gPad->SetLogy();
    }
    can_chani_calibcorr->Update();

    DisplayApp.Run();
}