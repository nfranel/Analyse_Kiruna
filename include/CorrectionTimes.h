#ifndef CORRECTIONTIMES_H
#define CORRECTIONTIMES_H

#include <cstdint>  // used for uint32_t
#include <Rtypes.h>
#include <TCanvas.h>

bool InWindow(const Double_t mid_interval, const Double_t half_window, const Double_t tested_value);
void file_exists(const std::string &filename);
std::vector<size_t> sort_index(const std::vector<Double_t>& vec);
void ADCFormating(const uint16_t *ener1, const uint16_t *ener2, Double_t &energy_adc, uint16_t *adc_channels, const std::string detector, const Double_t *correction, const uint16_t *pedestals);
Double_t ADCtoEnergyMaud(const Double_t ener_adc, const Double_t temp, const Double_t correction);
Double_t ADCtoEnergyDSSD(const uint16_t *adc_channels, Double_t *adc_channels_calib, const Double_t temp, const Double_t coeffs[64][8]);
void ADC_temp_correction(const uint16_t *adc_channels, Double_t *adc_channels_corr, const Double_t temp);
Double_t ADCtoEnergyCEA(const Double_t ener_adc, const Double_t temp);
void ExtractTemperatureValues(std::vector<Double_t> &tempvec, std::vector<uint32_t> &tsvec, const std::string timestamp_filename);
void FindTemperature(const uint32_t pps_timestamp, Double_t &temp_final, uint64_t &ite_init, const std::vector<uint32_t> timestamp_vec, std::vector<Double_t> temperature_vec);
void ExtractCorrectionValues(std::vector<Double_t> &corrvec, std::vector<uint32_t> &tsvec, const std::string correction_filename);
void FindCorrection(const uint32_t pps_timestamp, Double_t &correction, uint64_t &ite_init, const std::vector<uint32_t> timestamp_vec, std::vector<Double_t> correction_vec);
void ExtractDSSDCoefs(Double_t coeffs[64][8], const std::string dssd_p_coefs_filename, const std::string dssd_n_coefs_filename);
void FindPosition(uint16_t *adc_channels, Double_t *position, std::string detector);
void CombineDSSSDEvents(const std::string &init_tree_name, const std::string &final_tree_name);
void CorrectTimes(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &glitch_corr_count, uint32_t &minor_corr_count);
void ChangeTimeOriginIJCLAB(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &abs_gps_time_ref);
void ChangeTimeOriginUCD(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &abs_gps_time_ref);
void BuildFinalTreeIJCLAB(const std::string &init_tree_name, const std::string &final_tree_name, const std::vector<uint64_t> coinc_idx_ijclab, const std::vector<uint64_t> coinc_idx_cea, const std::vector<uint64_t> coinc_idx_ucda, const std::vector<uint64_t> coinc_idx_ucdb, const std::vector<uint64_t> coinc_idx_ucdc, const std::vector<uint64_t> coinc_idx_ucdd, std::string detector);
void BuildFinalTreeUCD(const std::string &tree_name_a, const std::string &tree_name_b, const std::string &tree_name_c, const std::string &tree_name_d, const std::string &final_tree_name, const std::vector<uint64_t> coinc_idx_ijclab1, const std::vector<uint64_t> coinc_idx_ijclab2, const std::vector<uint64_t> coinc_idx_cea, const std::vector<uint64_t> coinc_idx_a_b, const std::vector<uint64_t> coinc_idx_a_c, const std::vector<uint64_t> coinc_idx_a_d, const std::vector<uint64_t> coinc_idx_b_a, const std::vector<uint64_t> coinc_idx_b_c, const std::vector<uint64_t> coinc_idx_b_d, const std::vector<uint64_t> coinc_idx_c_a, const std::vector<uint64_t> coinc_idx_c_b, const std::vector<uint64_t> coinc_idx_c_d, const std::vector<uint64_t> coinc_idx_d_a, const std::vector<uint64_t> coinc_idx_d_b, const std::vector<uint64_t> coinc_idx_d_c);
void BuildFinalTreeCEA(const std::string &init_tree_name, const std::string &final_tree_name, const std::vector<uint64_t> coinc_idx_ijclab1, const std::vector<uint64_t> coinc_idx_ijclab2, const std::vector<uint64_t> coinc_idx_ucda, const std::vector<uint64_t> coinc_idx_ucdb, const std::vector<uint64_t> coinc_idx_ucdc, const std::vector<uint64_t> coinc_idx_ucdd);
void CombineUCDSubDets(const std::string &tree_name_a, const std::string &tree_name_b, const std::string &tree_name_c, const std::string &tree_name_d, const std::string &final_tree_name);
void FindCoincidences(const std::string &tree_name_a, const std::string &tree_name_b, std::vector<uint64_t> &idx_tree_a, std::vector<uint64_t> &idx_tree_b, std::vector<Double_t> &coinc_delays, const Double_t window_delay, const bool display, std::vector<TCanvas*> &save_canvas);
void DisplayDSSDChannels(const std::string &tree_name);

#endif