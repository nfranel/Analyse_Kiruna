#ifndef CORRECTIONTIMES_H
#define CORRECTIONTIMES_H

#include <cstdint>  // used for uint32_t
#include <Rtypes.h>

bool InWindow(const Double_t mid_interval, const Double_t half_window, const Double_t tested_value);
void file_exists(const std::string &filename);
void CorrectTimes(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &glitch_corr_count, uint32_t &minor_corr_count);
void ChangeTimeOrigin(const std::string &init_tree_name, const std::string &final_tree_name, std::string detector, uint32_t &abs_gps_time_ref);

#endif
