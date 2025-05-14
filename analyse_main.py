from Analyze_rootfiles import KirunaAnalysis
import matplotlib.pyplot as plt

# conda install -c conda-forge numpy matplotlib scipy pandas uproot memory_profiler pyqt

opening_mode = "compton"
# opening_mode = "spectro"
analysis = KirunaAnalysis(opening_mode=opening_mode)

# LC and spectra for all, ionisation maximum, burst and background
ergcut = None
detector = "all"
begins = ["begin", "begin", 56500, 130000]
ends = ["end", 15000, 74000, 320000]
xlog, ylog = False, True
for i in range(len(begins)):
    analysis.create_lightcurve(detector=detector, init_date=begins[i], end_date=ends[i], erg_cut=ergcut, bins=1000)
    analysis.create_spectrum(begins[i], ends[i], detector=detector, erg_cut=ergcut, xlog=xlog, ylog=ylog)
plt.show()


# Resolutions
# ~ 8% all
# ~ 9% ucd
# ~ 7% maud
#
# Selection en temps :
# total begin - end
# ionisation maximum begin - 15000
# burst 56500 - 74000
# bkg 130000 - 320000


# ergcut = None
# ergcut = (20, 1000)
# ergcut = (20, 200)
# periods = [["23/06/2024 4:20:00", "23/06/2024 18:10:00"],
#            ["24/06/2024 5:00:00", "24/06/2024 21:00:00"],
#            ["25/06/2024 6:40:00", "25/06/2024 22:50:00"],
#            ["26/06/2024 8:00:00", "end"]]
# periods = [["22/06/2024 23:20:00", "23/06/2024 00:20:00"]]
# init_date, end_date = "23/06/2024 11:50:00", "23/06/2024 15:26:40"
# init_date, end_date = "23/06/2024 11:50:00", "23/06/2024 11:55:00"
# analysis.create_spectrum(init_date, end_date)
# init_date, end_date = "23/06/2024 11:55:00", "23/06/2024 12:00:00"
# analysis.create_spectrum(init_date, end_date)
# init_date, end_date = "24/06/2024 00:00:00", "24/06/2024 00:05:00"
# analysis.create_spectrum(init_date, end_date)
# init_date, end_date = "24/06/2024 00:05:00", "24/06/2024 00:10:00"
# analysis.create_spectrum(init_date, end_date)

# analysis.fft_analysis(erg_cut=ergcut, mode="kiruna", add_signal=False, times_of_emptyness=False, periods=periods, figtitle="Day1")
# analysis.create_spectrum("begin", "end")
# analysis.create_spectrum("26/06/2024 8:00:00", "end")
# analysis.create_lightcurve(init_date="begin", end_date="end", erg_cut=ergcut)
# analysis.create_lightcurve(init_date="23/06/2024 11:50:00", end_date="23/06/2024 15:26:40", erg_cut=ergcut)


# save_ratio, save_date, save_r2, save_error, save_nev, save_reso = analysis.fit_511()
#
# with open("Kiruna_data/calib_maud/fine_calib.txt", "r") as f:
#     x_correc, y_correc = np.array([line.split(" ") for line in f.read().split("\n")[1:]], dtype=float).transpose()
#
# fit_time = (np.array([int(date.split(" ")[-1]) for date in save_date]) + np.array([int(date.split(" ")[0]) for date in save_date])) / 2
#
# nbins = 1000
# plt.rcParams.update({'font.size': 15})
# figglobal, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, constrained_layout=True)
#
# ax1.hist(analysis.d2b_df.time_corr_abs.values, bins=nbins, histtype="step", label=f"In flight light curve\nergcut = {None}")
# ax1.set(ylabel="Number of event")
#
# ax2.plot(x_correc, y_correc)
# ax2.set(ylabel="Calibration correction")
#
# ax3.plot(fit_time, save_reso)
# ax3.set(xlabel="Corrected time from 2024/06/22 20:00:00 (s)", ylabel="Fit resolution (%)")
#
# ax4.hist(save_reso, bins=30)
# ax4.set(xlabel="Fit resolution (%)", ylabel="Count")
#
# plt.show()
#
# print("val moy : ", np.mean(save_reso))


# fit_quad_calib(analysis.d2b_init_ts, analysis.d2b_final_ts, save_ratio, save_date, save_r2, save_error, save_nev)
# fit_quad_calib()

# timesbkg = ["24/06/2024 15:50:00", "24/06/2024 19:26:40"]
# timesburst = ["23/06/2024 11:50:00", "23/06/2024 15:26:40"]
# init1, end1 = analysis.get_timestamp(timesbkg[0]), analysis.get_timestamp(timesbkg[1])
# init2, end2 = analysis.get_timestamp(timesburst[0]), analysis.get_timestamp(timesburst[1])
#
# maud_df_selectbkg = analysis.d2b_df[np.logical_and(analysis.d2b_df.pps_cpt_corr_abs >= init1, analysis.d2b_df.pps_cpt_corr_abs <= end1)]
# maud_df_selectburst = analysis.d2b_df[np.logical_and(analysis.d2b_df.pps_cpt_corr_abs >= init2, analysis.d2b_df.pps_cpt_corr_abs <= end2)]
#
# bins = np.linspace(-100, 5500, 1000)
#
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
# count1, ener1 = ax1.hist(maud_df_selectbkg.energy.values, bins=bins, histtype="step", label="In flight spectrum bkg")[:2]
# count2, ener2 = ax1.hist(maud_df_selectburst.energy.values, bins=bins, histtype="step", label="In flight spectrum burst")[:2]
# ax1.axvline(511, label="511 keV", color="green")
# ax2.axvline(511, label="511 keV", color="green")
# ax2.plot((ener1[1:]+ener1[:-1])/2, count2 - count1, label="Burst - background")
# ax1.set(xlabel="Energy (keV)", ylabel="Number of event", yscale="log", ylim=(1, 300000))
# ax2.set(xlabel="Energy (keV)", ylabel="Number of event", yscale="log", xscale="log", ylim=(1, 300000))
#
# ax1.legend()
# ax2.legend()
# plt.show()


# # =================== Building the temperature file
# from datetime import datetime, timezone
#
# files = ["./Kiruna_data/Balloon_flight_complementary_data/CAEN_Temperature-data-2024-09-03 16_53_32.csv",
#          "./Kiruna_data/Balloon_flight_complementary_data/DSSD_Temperature-data-2024-09-03 16_53_15.csv",
#          "./Kiruna_data/Balloon_flight_complementary_data/MS5611_Temperature-data-2024-09-03 16_52_10.csv",
#          "./Kiruna_data/Balloon_flight_complementary_data/PI_Temperature-data-2024-09-03 16_52_49.csv"]
# savedfiles = ["./Kiruna_data/calib_maud/maud_time_vs_temp_flight.txt", "./Kiruna_data/calib_dssd/dssd_time_vs_temp_flight.txt",
#               "./Kiruna_data/Balloon_flight_complementary_data/MS5611_time_vs_temp_flight.txt",
#               "./Kiruna_data/Balloon_flight_complementary_data/PI_time_vs_temp_flight.txt"]
# for itefile in range(len(files)):
#     times = []
#     timestamps = []
#     timestamps_abs = []
#     temperatures = []
#     with open(files[itefile], "r") as f:
#         lines = f.read().split("\n")[2:]
#     for line in lines:
#         timeref = int(datetime(2024, 6, 22, 20, 0, 0).timestamp())
#         vals = line.split(",")
#         date, timeval = vals[0].split(" ")
#         year, month, day = map(int, date.split("-"))
#         hour, min, sec = map(int, timeval.split(":"))
#         ts = int(datetime(year, month, day, hour, min, sec).timestamp())
#         times.append(vals[0])
#         timestamps.append(ts)
#         timestamps_abs.append(ts - timeref)
#         temperatures.append(float(vals[1].split(" ")[0]))
#
#     with open(savedfiles[itefile], "w") as fs:
#         fs.write("LocalTime\tTimestamp\tTimestampNewRef\tTemperature\n")
#         for ite in range(len(times)):
#             fs.write(f"{times[ite]}\t{timestamps[ite]}\t{timestamps_abs[ite]}\t{temperatures[ite]}\n")




###########################################################################################################################################
#  DSSD calib analysis
###########################################################################################################################################

# calib = "./Kiruna_data/calib_dssd/nside_v2.calib"
# show_calib_params(calib)

# show_dssd_temp_dependance()
# show_dssd_temp_dependance(file_bulk="./Kiruna_data/calib_dssd/fit_corr_bulk.txt", file_peak="./Kiruna_data/calib_dssd/fit_corr_pedest.txt")