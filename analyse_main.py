# conda install -c conda-forge numpy matplotlib scipy pandas uproot memory_profiler pyqt


from Analyze_rootfiles import KirunaAnalysis, fit_maud_511, fit_quad_calib
import matplotlib.pyplot as plt
import numpy as np
# detlist = ["ucda", "maud"]
detlist = None
# opening_mode = "compton"
opening_mode = "full"
init_date = [2024, 6, 23, 12, 0, 0]
end_date = [2024, 6, 23, 16, 0, 0]
analysis = KirunaAnalysis(opening_mode=opening_mode, detlist=detlist, compressed=True, init_date=init_date, end_date=end_date)

# print(analysis.data.dtypes)


from Analyze_rootfiles import KirunaAnalysis, fit_maud_511, fit_quad_calib
import matplotlib.pyplot as plt
import numpy as np
detlist = ["ucda", "ucdb", "ucdc", "ucdd", "maud"]
# opening_mode = "compton"
opening_mode = "spectro"
analysis = KirunaAnalysis(opening_mode=opening_mode, detlist=detlist)

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

# #=================== Building the altitude file
# times = []
# timestamps = []
# timestamps_abs = []
# alts = []
# with open("./Kiruna_data/Balloon_flight_complementary_data/CAEN_Altitude_data-2024-09-03_16_53_32.csv", "r") as f:
#     lines = f.read().split("\n")[2:-1]
# for line in lines:
#     timeref = int(datetime(2024, 6, 22, 20, 0, 0).timestamp())
#     vals = line.split(";")
#     date, timeval = vals[0].split(" ")
#     day, month, year = map(int, date.split("/"))
#     hour, min, sec = map(int, timeval.split(":"))
#     ts = int(datetime(year, month, day, hour, min, sec).timestamp())
#     times.append(vals[0])
#     timestamps.append(ts)
#     timestamps_abs.append(ts - timeref)
#     alts.append(float(vals[1]))
#
# with open("./Kiruna_data/Balloon_flight_complementary_data/time_vs_alt_flight.txt", "w") as fs:
#     fs.write("LocalTime\tTimestamp\tTimestampNewRef\tAltitude\n")
#     for ite in range(len(times)):
#         fs.write(f"{times[ite]}\t{timestamps[ite]}\t{timestamps_abs[ite]}\t{alts[ite]}\n")




###########################################################################################################################################
#  DSSD calib analysis
###########################################################################################################################################
from Analyze_rootfiles import show_calib_params, show_dssd_temp_dependance
calib = "./Kiruna_data/calib_dssd/nside_v2.calib"
# show_calib_params(calib)

show_dssd_temp_dependance()
# show_dssd_temp_dependance(file_bulk="./Kiruna_data/calib_dssd/fit_corr_bulk.txt", file_peak="./Kiruna_data/calib_dssd/fit_corr_pedest.txt")



from Analyze_rootfiles import full_fit_maud_511
full_fit_maud_511(filename="corr_abs_maud_filevtemp.root", xlog="linear", ylog="linear")
full_fit_maud_511(filename="corr_abs_maud_filevtemp.root", xlog="linear", ylog="log")
full_fit_maud_511(filename="corr_abs_maud_file.root", xlog="linear", ylog="linear")
full_fit_maud_511(filename="corr_abs_maud_file.root", xlog="linear", ylog="log")

from Analyze_rootfiles import fit_maud_511, fit_quad_calib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
mpl.use("Qt5Agg")
# To do when the correction is set to 1 (not according to the file then)
# Obtain the ratio and all fitting results
# save_ratio1, save_date1, save_r21, save_error1, save_nev1, save_reso1, d2b_init_ts1, d2b_final_ts1 = fit_maud_511(filename="corr_abs_maud_filevuncorr.root", graphs=True, begin=None, end=19000)
# save_ratio2, save_date2, save_r22, save_error2, save_nev2, save_reso2, d2b_init_ts2, d2b_final_ts2 = fit_maud_511(filename="corr_abs_maud_filevuncorr.root", graphs=True, begin=19000, end=None)
# save_ratio, save_date = save_ratio1 + save_ratio2, save_date1 + save_date2
save_ratio, save_date, save_r2, save_error, save_nev, save_reso, d2b_init_ts, d2b_final_ts = fit_maud_511(filename="corr_abs_maud_filevuncorr.root", graphs=True)
fit_time = (np.array([int(date.split(" ")[-1]) for date in save_date]) + np.array([int(date.split(" ")[0]) for date in save_date])) / 2
# Saving the fitting results
with open("Kiruna_data/calib_maud/saved_calib_results_before_p1.txt", "w") as f:
    f.write("ratio|time|date|r2|error|nev|reso\n")
    for ite in range(len(save_ratio) - 1):
        f.write(f"{save_ratio[ite]}|{fit_time[ite]}|{save_date[ite]}|{save_r2[ite]}|{save_error[ite]}|{save_nev[ite]}|{save_reso[ite]}\n")
    f.write(f"{save_ratio[-1]}|{fit_time[-1]}|{save_date[ite]}|{save_r2[-1]}|{save_error[-1]}|{save_nev[-1]}|{save_reso[-1]}\n")

# Retrieving the data fitted and saved
save_ratio, save_date, save_r2, save_error, save_nev, save_reso, fit_time = [], [], [], [], [], [], []
with open("Kiruna_data/calib_maud/saved_calib_results_before.txt", "r") as f:
    data = f.read().split("\n")[1:-1]
for dat in data:
    val_ratio, val_time, val_date, val_r2, val_error, val_nev, val_reso = dat.split("|")
    save_ratio.append(float(val_ratio))
    fit_time.append(float(val_time))
    save_date.append(val_date)
    save_r2.append(float(val_r2))
    save_error.append(float(val_error))
    save_nev.append(float(val_nev))
    save_reso.append(float(val_reso))

val_init, val_final = 19265, 323625
smoothed_time, smoothed_ratio = fit_quad_calib(val_init, val_final, save_ratio, save_date, save_r2, save_error, save_nev)

ts_temp_list = []
temp_list = []
with open("Kiruna_data/calib_maud/maud_time_vs_temp_flight.txt", "r") as f:
    data = f.read().split("\n")[1:-1]
for dat in data:
    ts, temp = dat.split("\t")[2:]
    ts_temp_list.append(float(ts))
    temp_list.append(float(temp))

ts_alt_list = []
alt_list = []
with open("./Kiruna_data/Balloon_flight_complementary_data/time_vs_alt_flight.txt", "r") as f:
    data = f.read().split("\n")[1:-1]
for dat in data:
    ts, alt = dat.split("\t")[2:]
    ts_alt_list.append(float(ts))
    alt_list.append(float(alt))

fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
ax1.axhline(1, color="red")
ax1.plot(fit_time, save_ratio, color="blue", label=r"Temperature correction coefficient $C_{D2B}$")
ax1.legend()
ax1.set(ylabel=r"$C_{D2B}$")
ax2.plot(ts_temp_list, temp_list, color="blue", label="Temperature variation")
ax2.set(xlabel="Time (s)", ylabel="Temperature (°C)", xlim=(3600, None))
plt.tight_layout()
plt.show()

fig2, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
ax1.axhline(1, color="red")
ax1.plot(fit_time, save_ratio, color="blue", label=r"Temperature correction coefficient $C_{D2B}$")
ax1.legend()
ax1.set(ylabel=r"$C_{D2B}$")
ax2.plot(ts_temp_list, temp_list, color="blue", label="Temperature variation")
ax2.set(ylabel="Temperature (°C)", xlim=(3600, None))
ax3.plot(ts_alt_list, alt_list, color="blue", label="Altitude variation")
ax3.set(xlabel="Time (s)", ylabel="Altitude (m)", xlim=(3600, None))
plt.tight_layout()
plt.show()


# To do when the carroction was made, to verify if the correction is properly made.
save_ratio, save_date, save_r2, save_error, save_nev, save_reso, d2b_init_ts, d2b_final_ts = fit_maud_511(filename="corr_abs_maud_file.root", graphs=True, begin=19000, end=None)
fit_time = (np.array([int(date.split(" ")[-1]) for date in save_date]) + np.array([int(date.split(" ")[0]) for date in save_date])) / 2
with open("Kiruna_data/calib_maud/saved_calib_results_after.txt", "w") as f:
    f.write("ratio|time|date|r2|error|nev|reso\n")
    for ite in range(len(save_ratio) - 1):
        f.write(f"{save_ratio[ite]}|{fit_time[ite]}|{save_date[ite]}|{save_r2[ite]}|{save_error[ite]}|{save_nev[ite]}|{save_reso[ite]}\n")
    f.write(f"{save_ratio[-1]}|{fit_time[-1]}|{save_date[ite]}|{save_r2[-1]}|{save_error[-1]}|{save_nev[-1]}|{save_reso[-1]}\n")

fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))
ax1.axhline(1, color="red")
ax1.plot(fit_time, save_ratio, color="blue", label=r"Residual absolute error $511_{fit}$ / 511")
ax1.legend()
ax1.set(xlabel="Time (s)", ylabel=r"$511_{fit}$ / 511")
ax2.hist(save_reso, bins=30)
ax2.set(xlabel="Fit resolution (%)", ylabel="Number of fit performed")

plt.tight_layout()
plt.show()


#################################################################################################
# Light curves
#################################################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import uproot
from datetime import datetime, timezone

mpl.use("Qt5Agg")
plt.rcParams.update({'font.size': 13})
filenames = ["corr_abs_maud_file.root", "corr_abs_dssd_file.root", "corr_abs_ucda_file.root",
             "corr_abs_ucdb_file.root", "corr_abs_ucdc_file.root", "corr_abs_ucdd_file.root"]
detec = ["D2B", "D1B", "D2AA", "D2AB", "D2AC", "D2AD"]
enerlabels = []
lclabels = []

df_list = []
begin = 3600
end = "end"
for fileite, filename in enumerate(filenames):
    with uproot.open(filename) as file:
        # here, file is a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        saved_columns = ["time_corr_abs", "pps_cpt_corr_abs", "energy"]

        # convert the tree in dataframe
        data = pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=100000, library="pd")], ignore_index=True)

    init_time = np.min(data["time_corr_abs"].values)
    final_time = np.max(data["time_corr_abs"].values)

    if begin == "begin":
        begin_time = int(init_time)
    else:
        begin_time = begin
    if end == "end":
        end_time = int(final_time) + 1
    else:
        end_time = end

    df_list.append(data[(data.time_corr_abs >= begin_time) & (data.time_corr_abs <= end_time)])
    enerlabels.append(f"{detec[fileite]} spectrum over {begin_time}-{end_time} time range")
    lclabels.append(f"{detec[fileite]} light curve")


nbins = 1000
enerbins_scint = np.linspace(1, 5500, nbins)
enerbins_dssd = np.linspace(1, 1650, 300)
timebins = np.arange(begin_time, end_time+60, 60)
# Automatic all light curves
for ite in range(len(df_list)):
    if ite == 1:
        enerbins = enerbins_dssd
    else:
        enerbins = enerbins_scint
    # figener, ax = plt.subplots(1, 1, figsize=(10, 6))
    # ax.hist(df_list[ite].energy.values, bins=enerbins, histtype="step", label=enerlabels[ite])
    # ax.set(xlabel="Energy (keV)", ylabel="Number of events")
    # ax.legend()
    # plt.tight_layout()
    # plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/{detec[ite]}_spec_{begin_time}-{end_time}.png")
    # plt.close(figener)

    figenerlog, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.hist(df_list[ite].energy.values, bins=enerbins, histtype="step", label=enerlabels[ite])
    ax.set(xlabel="Energy (keV)", ylabel="Number of events", yscale="log")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/{detec[ite]}_speclog_{begin_time}-{end_time}.png")
    plt.close(figenerlog)


    def gauss_func(x, ampl, mu, sigma, a, b, x0):
        """

        :param x:
        :param ampl:
        :param mu:
        :param sigma:
        :param a:
        :param b:
        :param x0:
        :return:
        """
        return ampl * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + a * (x - x0) + b
    figenerlog, ax = plt.subplots(1, 1, figsize=(10, 6))
    h = ax.hist(df_list[ite].energy.values, bins=enerbins, histtype="step", label=enerlabels[ite])
    itei, itef = 80, 110
    centro = (h[1][itei:itef] + h[1][itei + 1:itef + 1]) / 2
    vals = h[0][itei:itef]
    popt, pcov = curve_fit(gauss_func, centro, vals, p0=[10000, 511, 30, 1, 1, 1])[:2]
    ax.plot(centro, gauss_func(centro, *popt))
    ax.set(xlabel="Energy (keV)", ylabel="Number of events", yscale="log")
    ax.axvline(511, color="green", label="511 keV")
    ax.legend()
    plt.tight_layout()
    # plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/{detec[ite]}_speclog_511_{begin_time}-{end_time}.png")
    # plt.close(figenerlog)

    figlc, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.hist(df_list[ite].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[ite])
    ax.set(xlabel="Time (s)", ylabel="count rate (couts/min)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/{detec[ite]}_LC_{begin_time}-{end_time}.png")
    plt.close(figlc)


    # figlclog, ax = plt.subplots(1, 1, figsize=(10, 6))
    # ax.hist(df_list[ite].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[ite])
    # ax.set(xlabel="Time (s)", ylabel="count rate (couts/min)", yscale="log")
    # ax.legend()
    # plt.tight_layout()
    # plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/{detec[ite]}_LClog_{begin_time}-{end_time}.png")
    # plt.close(figlclog)


fig3lc, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
ax1.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[0])
ax1.set(ylabel="count rate\n(couts/min)")
ax1.legend()

ax2.hist(df_list[1].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[1])
ax2.set(ylabel="count rate\n(couts/min)")
ax2.legend()

ax3.hist(np.concatenate((df_list[2].time_corr_abs.values, df_list[3].time_corr_abs.values, df_list[4].time_corr_abs.values, df_list[5].time_corr_abs.values)), bins=timebins, histtype="step", label=f"D2A light curve")
ax3.set(xlabel="Time (s)", ylabel="count rate\n(couts/min)")
ax3.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/triple_LC_{begin_time}-{end_time}.png")
plt.close(fig3lc)

ts_alt_list = []
alt_list = []
with open("./Kiruna_data/Balloon_flight_complementary_data/time_vs_alt_flight.txt", "r") as f:
    data = f.read().split("\n")[1:-1]
for dat in data:
    ts, alt = dat.split("\t")[2:]
    ts_alt_list.append(float(ts))
    alt_list.append(float(alt))

fig3lcalt, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 7), sharex=True)
ax1.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[0])
ax1.set(ylabel="count rate\n(couts/min)")
ax1.legend()

ax2.hist(df_list[1].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[1])
ax2.set(ylabel="count rate\n(couts/min)")
ax2.legend()

ax3.hist(np.concatenate((df_list[2].time_corr_abs.values, df_list[3].time_corr_abs.values, df_list[4].time_corr_abs.values, df_list[5].time_corr_abs.values)), bins=timebins, histtype="step", label=f"D2A light curve")
ax3.set(ylabel="count rate\n(couts/min)")
ax3.legend()

ax4.plot(ts_alt_list, alt_list, label="Altitude variation")
ax4.set(xlabel="Time (s)", ylabel="Altitude (m)", xlim=(0, None))
ax4.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/triple_LC_alt_{begin_time}-{end_time}.png")
# plt.close(fig3lcalt)

ts_alt_list = []
alt_list = []
with open("./Kiruna_data/Balloon_flight_complementary_data/time_vs_alt_flight.txt", "r") as f:
    data = f.read().split("\n")[1:-1]
for dat in data:
    ts, alt = dat.split("\t")[2:]
    ts_alt_list.append(float(ts))
    alt_list.append(float(alt))

fig2lcalt, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 5), sharex=True)
ax1.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[0])
ax1.set(ylabel="count rate\n(couts/min)")
ax1.legend()

ax2.hist(np.concatenate((df_list[2].time_corr_abs.values, df_list[3].time_corr_abs.values, df_list[4].time_corr_abs.values, df_list[5].time_corr_abs.values)), bins=timebins, histtype="step", label=f"D2A light curve")
ax2.set(ylabel="count rate\n(couts/min)")
ax2.legend()

ax3.plot(ts_alt_list, alt_list, label="Altitude variation")
ax3.set(xlabel="Time (s)", ylabel="Altitude (m)", xlim=(0, None))
ax3.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/double_LC_alt_{begin_time}-{end_time}.png")
# plt.close(fig2lcalt)



#################################################################################################
# Background
#################################################################################################
# Coordinates of the different position on the trajectory:
# lat1, lon1 = 67, 6
lat2, lon2 = 68, -16
# lat3, lon3 = 70, -36
# lat4, lon4 = 71, -57
# We keep   lat2, lon2 = 68, -16   as it is representative of the flight
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from datetime import datetime, timezone

mpl.use("Qt5Agg")
plt.rcParams.update({'font.size': 13})
filenames = ["corr_abs_maud_file.root", "corr_abs_dssd_file.root", "corr_abs_ucda_file.root", "corr_abs_ucdb_file.root", "corr_abs_ucdc_file.root", "corr_abs_ucdd_file.root"]
detec = ["D2B", "D1B", "D2AA", "D2AB", "D2AC", "D2AD"]
enerlabels = []
lclabels = []

df_list = []
begin = 155000
end = 245000
enermin = 400
enermax = 1000
for fileite, filename in enumerate(filenames):
    with uproot.open(filename) as file:
        # here, file is a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        saved_columns = ["time_corr_abs", "pps_cpt_corr_abs", "energy"]

        # convert the tree in dataframe
        data = pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=100000, library="pd")], ignore_index=True)

    init_time = np.min(data["time_corr_abs"].values)
    final_time = np.max(data["time_corr_abs"].values)

    if begin == "begin":
        begin_time = int(init_time)
    else:
        begin_time = begin
    if end == "end":
        end_time = int(final_time) + 1
    else:
        end_time = end

    temp_df = data[(data.time_corr_abs >= begin_time) & (data.time_corr_abs <= end_time)]
    df_list.append(temp_df[(temp_df.energy >= enermin) & (temp_df.energy <= enermax)])
    enerlabels.append(f"{detec[fileite]} spectrum over {begin_time}-{end_time} time range")
    lclabels.append(f"{detec[fileite]} light curve")

# Results :
print("=====================================================================================")
print(f" In the {enermin} keV - {enermax} keV energy range we obtain between {begin} s and {end} s :")
print(f"      {len(df_list[1])} events, making {np.around(len(df_list[1]) / (end - begin), 2)} counts/s in D1B")
print(f"      {len(df_list[0])} events, making {np.around(len(df_list[0]) / (end - begin), 2)} counts/s in D2B")
print(f"    {len(df_list[2])} events, making {np.around(len(df_list[2]) / (end - begin), 2)} counts/s in D2AA")
print(f"    {len(df_list[3])} events, making {np.around(len(df_list[3]) / (end - begin), 2)} counts/s in D2AB")
print(f"    {len(df_list[4])} events, making {np.around(len(df_list[4]) / (end - begin), 2)} counts/s in D2AC")
print(f"    {len(df_list[5])} events, making {np.around(len(df_list[5]) / (end - begin), 2)} counts/s in D2AD")
print(f"      {len(df_list[2]) + len(df_list[3]) + len(df_list[4]) + len(df_list[5])} events, making {np.around((len(df_list[2]) + len(df_list[3]) + len(df_list[4]) + len(df_list[5])) / (end - begin), 2)} counts/s in D2A")
print(f"         {len(df_list[0]) + len(df_list[2]) + len(df_list[3]) + len(df_list[4]) + len(df_list[5])} events, making {np.around((len(df_list[0]) + len(df_list[2]) + len(df_list[3]) + len(df_list[4]) + len(df_list[5])) / (end - begin), 2)} counts/s in D2A and D2B")
print(f"         {len(df_list[0]) + len(df_list[1]) + len(df_list[2]) + len(df_list[3]) + len(df_list[4]), + len(df_list[5])} events, making {np.around((len(df_list[0]) + len(df_list[1]) + len(df_list[2]) + len(df_list[3]) + len(df_list[4]) + len(df_list[5])) / (end - begin), 2)} counts/s in D1B, D2A and D2B")


# Analysis of simulated background
# to be done in the terminal with megalib activated and
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
import sys
# sys.path.append("/home/nathan/Desktop/GRB_polarimetry/Simu/GRB-simulation/Sim_code/grb-COMCUBE-simu")
from src.General.funcmod import analyze_bkg_event, det_counter

def gauss(x, amp, mu, sig):
  """

  """
  return amp * norm.pdf(x, loc=mu, scale=sig)

ergcut = (30, 1000)
data_prefix = "/home/nathan/Desktop/GRB_polarimetry/Simu/GRB-simulation/Sim_code/COMCUBEBAL/rawsim/"
data_files = [f"{data_prefix}baloonbkg_Background-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Alpha-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Electron-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Muonm-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Muonp-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Neutron-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Photon-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Positron-sat0-0000_0.0_0.0.inc1.id1.extracted.tra",
              f"{data_prefix}baloonbkg_Background_Proton-sat0-0000_0.0_0.0.inc1.id1.extracted.tra"]
lat = 70
alt = 40
geometry = "/home/nathan/Desktop/GRB_polarimetry/Simu/GRB-simulation/Sim_code/COMCUBEBAL/v2COMCUBE-ballon-2023.geo.setup"
array_dtype = np.float64
bkgsinglecont = []
particle = ["Full background", "Alpha particles", "Electrons", "Muons-", "Muons+", "Neutrons", "Photons", "Positrons", "Protons"]
print("===========================================================")
print(f"  Simulated detector count rate   -   ergcut : {ergcut}")
print("===========================================================")
for data_ite, data_file in enumerate(data_files):
    decbkg, altbkg, compton_second, compton_ener, compton_time, single_ener, single_time, compton_first_detector, compton_sec_detector, single_detector = analyze_bkg_event(data_file, lat, alt, geometry, array_dtype)
    df_compton = pd.DataFrame({"compton_ener": compton_ener, "compton_second": compton_second, "compton_time": compton_time,
                               "compton_first_detector": compton_first_detector, "compton_sec_detector": compton_sec_detector})
    df_single = pd.DataFrame({"single_ener": single_ener, "single_time": single_time, "single_detector": single_detector})
    df_compton = df_compton[(df_compton.compton_ener >= ergcut[0]) & (df_compton.compton_ener <= ergcut[1])]
    df_single = df_single[(df_single.single_ener >= ergcut[0]) & (df_single.single_ener <= ergcut[1])]
    bkgsinglecont.append(df_single)
    det_stat_compton = det_counter(np.concatenate((df_compton.compton_first_detector.values, df_compton.compton_sec_detector.values))).flatten()
    det_stat_single = det_counter(df_single.single_detector.values).flatten()

    d1a_cr = (det_stat_compton[3] + det_stat_single[3]) / 3600
    d1b_cr = (det_stat_compton[2] + det_stat_single[2]) / 3600
    d2a_cr = (det_stat_compton[4] + det_stat_single[4]) / 3600
    d2b_cr = (det_stat_compton[0] + det_stat_single[0]) / 3600

    print(f"   = {particle[data_ite]} = ")
    print(f"     D1A : {d1a_cr} hit/s")
    print(f"     D1B : {d1b_cr} hit/s")
    print(f"     D2A : {d2a_cr} hit/s")
    print(f"     D2B : {d2b_cr} hit/s")


df_selec = bkgsinglecont[0][bkgsinglecont[0].single_detector == 3]

plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
h = np.histogram(df_selec.single_ener, bins=100)
try:
    cent = h[1][45:91]
    vals = h[0][45:90]
    centroids = (cent[1:] + cent[:-1]) / 2
    # plt.axvline(cent[0])
    # plt.axvline(cent[-1])
    popt, pcov = curve_fit(gauss, centroids, np.log10(vals), bounds=((-np.inf, 500, 0), (np.inf, 900, np.inf)))[:2]
    ax.plot(centroids, 10**gauss(centroids, *popt), label="Gaussian fit of the bump")
except:
    pass
ax.stairs(h[0], h[1])
ax.set(yscale="log", xlabel="Energy (keV)", ylabel="Number of events")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/bkg/fullbkg_spec_D1B_{ergcut[0]}-{ergcut[1]}")
plt.close(fig)
print(popt)


for data_ite in range(1, len(particle)):
    df_selec = bkgsinglecont[data_ite][bkgsinglecont[data_ite].single_detector == 3]
    plt.rcParams.update({'font.size': 15})
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    h = np.histogram(df_selec.single_ener, bins=100)
    if data_ite in [2, 8]:
        try:
            cent = h[1][45:91]
            vals = h[0][45:90]
            centroids = (cent[1:] + cent[:-1]) / 2
            # plt.axvline(cent[0])
            # plt.axvline(cent[-1])
            popt, pcov = curve_fit(gauss, centroids, np.log10(vals), bounds=((-np.inf, 500, 0), (np.inf, 900, np.inf)))[:2]
            print(popt)
            ax.plot(centroids, 10 ** gauss(centroids, *popt), label="Gaussian fit of the bump")
        except:
            pass
    ax.stairs(h[0], h[1], label=f"{particle[data_ite]} background spectrum")
    ax.set(yscale="log", xlabel="Energy (keV)", ylabel="Number of events")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/bkg/{particle[data_ite]}_spec_D1B_{ergcut[0]}-{ergcut[1]}")
    plt.close(fig)


df_selec = bkgsinglecont[0][bkgsinglecont[0].single_detector == 1]

plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
h = np.histogram(df_selec.single_ener, bins=100)
# cent = h[1][45:91]
# vals = h[0][45:90]
# centroids = (cent[1:] + cent[:-1]) / 2
# plt.axvline(cent[0])
# plt.axvline(cent[-1])
# popt, pcov = curve_fit(gauss, centroids, np.log10(vals), bounds=((-np.inf, 500, 0), (np.inf, 900, np.inf)))[:2]
# ax.plot(centroids, 10**gauss(centroids, *popt), label="Gaussian fit of the bump")
ax.stairs(h[0], h[1])
ax.set(yscale="log", xlabel="Energy (keV)", ylabel="Number of events")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/bkg/fullbkg_spec_D2B_{ergcut[0]}-{ergcut[1]}")
plt.close(fig)

df_selec = bkgsinglecont[0][bkgsinglecont[0].single_detector == 5]

plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
h = np.histogram(df_selec.single_ener, bins=100)
# cent = h[1][45:91]
# vals = h[0][45:90]
# centroids = (cent[1:] + cent[:-1]) / 2
# plt.axvline(cent[0])
# plt.axvline(cent[-1])
# popt, pcov = curve_fit(gauss, centroids, np.log10(vals), bounds=((-np.inf, 500, 0), (np.inf, 900, np.inf)))[:2]
# ax.plot(centroids, 10**gauss(centroids, *popt), label="Gaussian fit of the bump")
ax.stairs(h[0], h[1])
ax.set(yscale="log", xlabel="Energy (keV)", ylabel="Number of events")
# ax.axvline(511, color="green", label="511 keV")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/bkg/fullbkg_spec_D2A_{ergcut[0]}-{ergcut[1]}")
plt.close(fig)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from datetime import datetime, timezone

mpl.use("Qt5Agg")
plt.rcParams.update({'font.size': 15})
filenames = ["corr_abs_dssd_file.root"]
detec = ["D1B"]
enerlabels = []
lclabels = []

df_list = []
begin = 3600
end = "end"
for fileite, filename in enumerate(filenames):
    with uproot.open(filename) as file:
        # here, file is a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        saved_columns = ["time_corr_abs", "energy", "temperature"]

        # convert the tree in dataframe
        data = pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=100000, library="pd")], ignore_index=True)

    init_time = np.min(data["time_corr_abs"].values)
    final_time = np.max(data["time_corr_abs"].values)

    if begin == "begin":
        begin_time = int(init_time)
    else:
        begin_time = begin
    if end == "end":
        end_time = int(final_time) + 1
    else:
        end_time = end

    df_list.append(data[(data.time_corr_abs >= begin_time) & (data.time_corr_abs <= end_time)])
    enerlabels.append(f"{detec[fileite]} spectrum over {begin_time}-{end_time} time range")
    lclabels.append(f"{detec[fileite]} light curve")

# print centre of bins :
dssd_df = df_list[0]
print(np.around(np.mean(dssd_df[(dssd_df.temperature > 10)].temperature.values), 2))
print(np.around(np.mean(dssd_df[(dssd_df.temperature <= 10) & (dssd_df.temperature > 5)].temperature.values), 2))
print(np.around(np.mean(dssd_df[(dssd_df.temperature <= 5) & (dssd_df.temperature > 0)].temperature.values), 2))
print(np.around(np.mean(dssd_df[(dssd_df.temperature <= 0) & (dssd_df.temperature > -5)].temperature.values), 2))
print(np.around(np.mean(dssd_df[(dssd_df.temperature <= -5) & (dssd_df.temperature > -10)].temperature.values), 2))
print(np.around(np.mean(dssd_df[(dssd_df.temperature <= -10)].temperature.values), 2))

def calib_calculator(e, coefa, coefb, coefc, coefd, x1, x2, x3, x4, t1, t2, t3, t4):
    # print(e, coefa, coefb, coefc, coefd, x1, x2, x3, x4, t1, t2, t3, t4)
    f1 = (e - coefa - coefb * x1 - coefc * x1**2 - coefd * x1**3) / t1
    f2 = (e - coefa - coefb * x2 - coefc * x2**2 - coefd * x2**3) / t2
    f3 = (e - coefa - coefb * x3 - coefc * x3**2 - coefd * x3**3) / t3
    f4 = (e - coefa - coefb * x4 - coefc * x4**2 - coefd * x4**3) / t4
    # print(f1)
    # print(f2)
    # print(f3)
    # print(f4)

    f12 = (f2 - f1) / (x2 - x1)
    f13 = (f3 - f1) / (x3 - x1)
    f14 = (f4 - f1) / (x4 - x1)
    xca12 = (x2**2 - x1**2) / (x2 - x1)
    xca13 = (x3**2 - x1**2) / (x3 - x1)
    xca14 = (x4**2 - x1**2) / (x4 - x1)
    xcu12 = (x2**3 - x1**3) / (x2 - x1)
    xcu13 = (x3**3 - x1**3) / (x3 - x1)
    xcu14 = (x4**3 - x1**3) / (x4 - x1)

    f123 = (f13 - f12) / (xca13 - xca12)
    f124 = (f14 - f12) / (xca14 - xca12)
    x123 = (xcu13 - xcu12) / (xca13 - xca12)
    x124 = (xcu14 - xcu12) / (xca14 - xca12)

    d = (f124 - f123) / (x124 - x123)
    c = f123 - d * x123
    b = f12 - c * xca12 - d * xcu12
    a = f1 - b * x1 - c * x1**2 - d * x1**3

    if ~((np.around(f1, 5) == np.around(a + b * x1 + c * x1**2 + d * x1**3, 5)) &
          (np.around(f2, 5) == np.around(a + b * x2 + c * x2**2 + d * x2**3, 5)) &
          (np.around(f3, 5) == np.around(a + b * x3 + c * x3**2 + d * x3**3, 5)) &
          (np.around(f4, 5) == np.around(a + b * x4 + c * x4**2 + d * x4**3, 5))):
        raise ValueError
    print(a, b, c, d)
    return a, b, c, d

with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/nside_v3.calib") as f:
    coefs25 = f.read().split("\n")[:-1]
with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_bump.txt") as f:
    bumps_chan = f.read().split("\n")[:-1]
temps = np.array([14.09, 7.26, 2.62, -2.55, -7.56, -13.14]) - 25
e = 704 / 1000
combinations = [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5), (0, 1, 3, 4), (0, 1, 3, 5), (0, 1, 4, 5), (0, 2, 3, 4), (0, 2, 3, 5), (0, 2, 4, 5), (0, 3, 4, 5), (1, 2, 3, 4), (1, 2, 3, 5), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)]
coefs = np.zeros((32, 4))
coefs_var = np.zeros((32, 4, 15))
for i in range(32):
    dats = coefs25[i].split("\t")
    bumps = np.array(bumps_chan[i].split("|")[0].split(" "), dtype=float) / 1000
    coefa, coefb, coefc, coefd = float(dats[1]), float(dats[3]), float(dats[5]), float(dats[7])
    coefs[i] = np.array([coefa, coefb, coefc, coefd])
    for itecomb, comb in enumerate(combinations):
        print(coefa, coefb, coefc, coefd)
        vara, varb, varc, vard = calib_calculator(e, coefa, coefb, coefc, coefd, float(bumps[comb[0]]), float(bumps[comb[1]]),
                                                  float(bumps[comb[2]]),float(bumps[comb[3]]), temps[comb[0]], temps[comb[1]],
                                                  temps[comb[2]], temps[comb[3]])
        coefs_var[i][0][itecomb] = vara
        coefs_var[i][1][itecomb] = varb
        coefs_var[i][2][itecomb] = varc
        coefs_var[i][3][itecomb] = vard

    eneradc = bumps[0]
    tempener = temps[0]
    atest = coefs[i][0] + np.mean(coefs_var[i][0]) * tempener
    btest = coefs[i][1] + np.mean(coefs_var[i][1]) * tempener
    ctest = coefs[i][2] + np.mean(coefs_var[i][2]) * tempener
    dtest = coefs[i][3] + np.mean(coefs_var[i][3]) * tempener

    enercalc = 1000 * (atest + btest * eneradc + ctest * eneradc**2 + dtest * eneradc**3)
    print(enercalc)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

mpl.use("Qt5Agg")
def affine(x, a, b):
    return a*x+b

temperatures = np.array([14.09, 7.26, 2.62, -2.55, -7.56, -13.14])

with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_raw_bump.txt") as f:
    data = f.read().split("\n")[:-1]
    raw_bumps = [np.array(dat.split("|")[0].split(" "), dtype=float) for dat in data]
    raw_bumps_err = [np.array(dat.split("|")[1].split(" "), dtype=float) for dat in data]
rawfit_a = []
rawfit_b = []
fig_raw, axes1 = plt.subplots(4, 8, figsize=(25, 15))
for ite in range(32):
    popt, pcov = curve_fit(affine, temperatures, raw_bumps[ite], sigma=raw_bumps_err[ite])[:2]
    rawfit_a.append(popt[0])
    rawfit_b.append(popt[1])
    axes1[int(ite / 8)][ite % 8].errorbar(temperatures, raw_bumps[ite], yerr=raw_bumps_err[ite], color="blue")
    axes1[int(ite / 8)][ite % 8].plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red")
    axes1[int(ite / 8)][ite % 8].set(title=f"{round(popt[0], 2)} x + {round(popt[1], 2)}", xlabel="Temperature",
                                     ylabel="Peak position(ADC)")
fig_raw.set_constrained_layout(True)
plt.savefig("/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/dssd_calibration/bump_variation_fits_chanraw_all.png")
plt.close(fig_raw)
# plt.show()

with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_calib_bump.txt") as f:
    data = f.read().split("\n")[:-1]
    calib_bumps = [np.array(dat.split("|")[0].split(" "), dtype=float) for dat in data]
    calib_bumps_err = [np.array(dat.split("|")[1].split(" "), dtype=float) for dat in data]
calibfit_a = []
calibfit_b = []
fig_calib, axes1 = plt.subplots(4, 8, figsize=(25, 15))
for ite in range(32):
    popt, pcov = curve_fit(affine, temperatures, calib_bumps[ite], sigma=calib_bumps_err[ite])[:2]
    calibfit_a.append(popt[0])
    calibfit_b.append(popt[1])
    axes1[int(ite / 8)][ite % 8].errorbar(temperatures, calib_bumps[ite], yerr=calib_bumps_err[ite], color="blue")
    axes1[int(ite / 8)][ite % 8].plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red")
    axes1[int(ite / 8)][ite % 8].set(title=f"{round(popt[0], 2)} x + {round(popt[1], 2)}", xlabel="Temperature",
                                     ylabel="Peak position(ADC)")
fig_calib.set_constrained_layout(True)
plt.savefig("/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/dssd_calibration/bump_variation_fits_chancalib_all.png")
plt.close(fig_calib)

with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_calib_corr_bump.txt") as f:
    data = f.read().split("\n")[:-1]
    calib_corr_bumps = [np.array(dat.split("|")[0].split(" "), dtype=float) for dat in data]
    calib_corr_bumps_err = [np.array(dat.split("|")[1].split(" "), dtype=float) for dat in data]
calibcorrfit_a = []
calibcorrfit_b = []
fig_calib_corr, axes1 = plt.subplots(4, 8, figsize=(25, 15))
for ite in range(32):
    popt, pcov = curve_fit(affine, temperatures, calib_corr_bumps[ite], sigma=calib_corr_bumps_err[ite])[:2]
    calibcorrfit_a.append(popt[0])
    calibcorrfit_b.append(popt[1])
    axes1[int(ite / 8)][ite % 8].errorbar(temperatures, calib_corr_bumps[ite], yerr=calib_corr_bumps_err[ite], color="blue")
    axes1[int(ite / 8)][ite % 8].plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red")
    axes1[int(ite / 8)][ite % 8].set(title=f"{round(popt[0], 2)} x + {round(popt[1], 2)}", xlabel="Temperature",
                                     ylabel="Peak position(ADC)")
fig_calib_corr.set_constrained_layout(True)
plt.savefig("/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/dssd_calibration/bump_variation_fits_chancalibcorr_all.png")
plt.close(fig_calib_corr)



import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

mpl.use("Qt5Agg")
def affine(x, a, b):
    return a*x+b

temperatures = np.array([14.09, 7.26, 2.62, -2.55, -7.56, -13.14])
plt.rcParams.update({'font.size': 15})
ite_chan = 20

with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_raw_bump.txt") as f:
    data = f.read().split("\n")[:-1]
    raw_bumps = [np.array(dat.split("|")[0].split(" "), dtype=float) for dat in data]
    raw_bumps_err = [np.array(dat.split("|")[1].split(" "), dtype=float) for dat in data]
fig_fitchanraw, ax1 = plt.subplots(1, 1, figsize=(10, 6))
popt, pcov = curve_fit(affine, temperatures, raw_bumps[ite_chan], sigma=raw_bumps_err[ite_chan])[:2]
ax1.errorbar(temperatures, raw_bumps[ite_chan], yerr=raw_bumps_err[ite_chan], color="blue")
ax1.plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red",
         label=f"Fitted variation of the bump : {round(popt[0], 2)} x + {round(popt[1], 2)}")
ax1.legend()
ax1.set(xlabel="Temperature °C", ylabel="Bump position (ADC)")
fig_fitchanraw.set_constrained_layout(True)
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/dssd_calibration/bump_variation_fits_chan{ite_chan}raw_all.png")
plt.close(fig_fitchanraw)


with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_calib_bump.txt") as f:
    data = f.read().split("\n")[:-1]
    calib_bumps = [np.array(dat.split("|")[0].split(" "), dtype=float) for dat in data]
    calib_bumps_err = [np.array(dat.split("|")[1].split(" "), dtype=float) for dat in data]
fig_fitchancalib, ax2 = plt.subplots(1, 1, figsize=(10, 6))
popt, pcov = curve_fit(affine, temperatures, calib_bumps[ite_chan], sigma=calib_bumps_err[ite_chan])[:2]
ax2.errorbar(temperatures, calib_bumps[ite_chan], yerr=calib_bumps_err[ite_chan], color="blue")
ax2.plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red",
         label=f"Fitted variation of the bump : {round(popt[0], 2)} x + {round(popt[1], 2)}")
ax2.legend()
ax2.set(xlabel="Temperature °C", ylabel="Bump position (ADC)")
fig_fitchancalib.set_constrained_layout(True)
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/dssd_calibration/bump_variation_fits_chan{ite_chan}calib_all.png")
plt.close(fig_fitchancalib)

with open("/home/nathan/Desktop/GRB_polarimetry/Simu/Kiruna/Analyse_Kiruna/Kiruna_data/calib_dssd/fit_calib_corr_bump.txt") as f:
    data = f.read().split("\n")[:-1]
    calib_corr_bumps = [np.array(dat.split("|")[0].split(" "), dtype=float) for dat in data]
    calib_corr_bumps_err = [np.array(dat.split("|")[1].split(" "), dtype=float) for dat in data]
fig_fitchancalibcorr, ax3 = plt.subplots(1, 1, figsize=(10, 6))
popt, pcov = curve_fit(affine, temperatures, calib_corr_bumps[ite_chan], sigma=calib_corr_bumps_err[ite_chan])[:2]
ax3.errorbar(temperatures, calib_corr_bumps[ite_chan], yerr=calib_corr_bumps_err[ite_chan], color="blue")
ax3.plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red",
         label=f"Fitted variation of the bump : {round(popt[0], 2)} x + {round(popt[1], 2)}")
ax3.legend()
ax3.set(xlabel="Temperature °C", ylabel="Bump position (ADC)")
fig_fitchancalibcorr.set_constrained_layout(True)
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/dssd_calibration/bump_variation_fits_chan{ite_chan}calibcorr_all.png")
plt.close(fig_fitchancalibcorr)



#################################################################################################
# Light curves
#################################################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import uproot
from datetime import datetime, timezone

def date_to_exact_frame(date_info):
    timeref = int(datetime(2024, 6, 22, 18, 0, 0).timestamp())
    date, timeval = date_info.split(" ")
    day, month, year = map(int, date.split("/"))
    hour, min, sec = map(int, timeval.split(":"))
    ts = int(datetime(year, month, day, hour, min, sec).timestamp())
    return ts - timeref

# List of potential bursts :
grbdate = ["25/06/2024 22:58:35", "25/06/2024 7:32:40", "24/06/2024 19:19:45", "24/06/2024 6:01:11"]
grbtime = [1719356315, 1719300760, 1719256785, 1719208871]
grbtime_abs = [277115, 221560, 177585, 129671]
grbra = [90, 119, 240, 107]
grbdec = [-11, -5, -7, -9]
gondola_lon = [-58.9, -40.90, -30.74, -14.03]
gondola_lat = [71.65, 70.29, 70.15, 68.90]
grbvisible = ["~Nope", "Nope", "~Nope", "Nope"]
def radecconv(ra, dec):
    ras = np.array(ra.split(" "), dtype=float)
    decs = np.array(dec.split(" "), dtype=float)
    if decs[0] < 0:
        mult = -1
    else:
        mult = 1
    return (ras[0] + ras[1]/60 + ras[2]/3600)*15, (abs(decs[0]) + decs[1]/60 + decs[2]/3600) * mult
print(radecconv("06 00 04.8","-11 07 12"))
print(radecconv("07 54 57.6","-04 56 24"))
print(radecconv("15 53 50.4","-06 23 24"))
print(radecconv("06 49 46.6","-06 32 35"))

mpl.use("Qt5Agg")
plt.rcParams.update({'font.size': 13})
filenames = ["corr_abs_maud_file.root", "corr_abs_dssd_file.root", "corr_abs_ucda_file.root",
             "corr_abs_ucdb_file.root", "corr_abs_ucdc_file.root", "corr_abs_ucdd_file.root"]
detec = ["D2B", "D1B", "D2AA", "D2AB", "D2AC", "D2AD"]
enerlabels = []
lclabels = []

df_list = []
begin = 3600
end = "end"
for fileite, filename in enumerate(filenames):
    with uproot.open(filename) as file:
        # here, file is a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        saved_columns = ["time_corr_abs", "pps_cpt_corr_abs", "energy"]

        # convert the tree in dataframe
        data = pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=100000, library="pd")], ignore_index=True)

    init_time = np.min(data["time_corr_abs"].values)
    final_time = np.max(data["time_corr_abs"].values)

    if begin == "begin":
        begin_time = int(init_time)
    else:
        begin_time = begin
    if end == "end":
        end_time = int(final_time) + 1
    else:
        end_time = end

    df_list.append(data[(data.time_corr_abs >= begin_time) & (data.time_corr_abs <= end_time)])
    enerlabels.append(f"{detec[fileite]} spectrum over {begin_time}-{end_time} time range")
    lclabels.append(f"{detec[fileite]} light curve")

for ite_grb in range(len(grbtime_abs)):
    init_grbtime, end_grb_time = grbtime_abs[ite_grb] - 100, grbtime_abs[ite_grb] + 100
    timebins = np.arange(init_grbtime, end_grb_time+1, 1)
    # Automatic all light curves
    figlcgrb, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True)
    ax1.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[0])
    ax1.set(ylabel="count rate\n(couts/s)")
    ax1.legend()

    ax2.hist(df_list[1].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[1])
    ax2.set(ylabel="count rate\n(couts/s)")
    ax2.legend()

    ax3.hist(np.concatenate((df_list[2].time_corr_abs.values, df_list[3].time_corr_abs.values, df_list[4].time_corr_abs.values, df_list[5].time_corr_abs.values)), bins=timebins, histtype="step", label=f"D2A light curve")
    ax3.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
    ax3.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/triple_LC_GRB{ite_grb+1}_{init_grbtime}-{end_grb_time}.png")
    plt.close(figlcgrb)

plt.rcParams.update({'font.size': 15})
grblabel = ["Light curve around trigger time\nof GRB 240625957 ", "Light curve around trigger time\nof GRB 240625314 ", "Light curve around trigger time\nof GRB 240624805 ", "Light curve around trigger time\nof GRB 240624251"]
timebins1 = np.arange(grbtime_abs[0] - 100, grbtime_abs[0] + 100, 1)
timebins2 = np.arange(grbtime_abs[1] - 100, grbtime_abs[1] + 100, 1)
timebins3 = np.arange(grbtime_abs[2] - 100, grbtime_abs[2] + 100, 1)
timebins4 = np.arange(grbtime_abs[3] - 100, grbtime_abs[3] + 100, 1)
figlcgrb, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 8))
ax1.hist(df_list[0].time_corr_abs.values, bins=timebins1, histtype="step", label=grblabel[ite_grb])
ax1.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax1.legend()

ax2.hist(df_list[0].time_corr_abs.values, bins=timebins2, histtype="step", label=grblabel[ite_grb])
ax2.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax2.legend()

ax3.hist(df_list[0].time_corr_abs.values, bins=timebins3, histtype="step", label=grblabel[ite_grb])
ax3.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax3.legend()

ax4.hist(df_list[0].time_corr_abs.values, bins=timebins4, histtype="step", label=grblabel[ite_grb])
ax4.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax4.legend()

plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/D2B_LC_GRBs.png")
plt.close(figlcgrb)

# GOES X-ray satellite Solar flares
sf_dates_init = ["22/06/2024 19:34:00",
                 "23/06/2024 06:16:00", "23/06/2024 11:26:00", "23/06/2024 12:51:00",
                 "24/06/2024 04:08:00", "24/06/2024 04:45:00", "24/06/2024 11:09:00", "24/06/2024 11:43:00", "24/06/2024 19:02:00",
                 "25/06/2024 12:26:00"]
sf_dates_end = ["22/06/2024 20:04:00",
                "23/06/2024 06:48:00", "23/06/2024 11:50:00", "23/06/2024 13:11:00",
                "24/06/2024 04:23:00", "24/06/2024 04:56:00", "24/06/2024 11:15:00", "24/06/2024 11:53:00", "24/06/2024 19:14:00",
                "25/06/2024 13:06:00"]
sf_times_init_goes = [5640, 44160, 62760, 67860, 122880, 125100, 148140, 150180, 176520, 239160]
sf_times_end_goes = [7440, 46080, 64200, 69060, 123780, 125760, 148500, 150780, 177240, 241560]
sf_names = ["GOES_SF1", "GOES_SF2", "GOES_SF3", "GOES_SF4", "GOES_SF5",
            "GOES_SF6", "GOES_SF7", "GOES_SF8", "GOES_SF9", "GOES_SF10"]
for ite_sf in range(len(sf_times_init_goes)):
    init_sftime, end_sf_time = sf_times_init_goes[ite_sf] - 10 * 60, sf_times_end_goes[ite_sf] + 10 * 60
    timebins = np.arange(init_sftime, end_sf_time+1, 1)
    # Automatic all light curves
    figlcsf, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True)
    h1 = ax1.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[0])
    ax1.axvline(sf_times_init_goes[ite_sf])
    ax1.axvline(sf_times_end_goes[ite_sf])
    ax1.set(ylabel="count rate\n(couts/s)")
    ax1.legend()

    ax2.hist(df_list[1].time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[1])
    h2 = ax2.set(ylabel="count rate\n(couts/s)")
    ax2.legend()

    h3 = ax3.hist(np.concatenate((df_list[2].time_corr_abs.values, df_list[3].time_corr_abs.values, df_list[4].time_corr_abs.values, df_list[5].time_corr_abs.values)), bins=timebins, histtype="step", label=f"D2A light curve")
    ax3.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
    ax3.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/triple_LC_solarflare{ite_sf+1}_{init_sftime}-{end_sf_time}.png")
    plt.close(figlcsf)


timebinstot = np.arange(begin_time, end_time+60, 60)
colors = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan"]

plt.rcParams.update({'font.size': 13})
figfullsf = plt.figure(figsize=(10, 11))
gs = gridspec.GridSpec(6, 2, figure=figfullsf)  # 6 lignes, 2 colonnes

# First line, figure of the whole light curve
ax_big = figfullsf.add_subplot(gs[0, :])  # toute la 1ère ligne
ax_big.hist(df_list[0].time_corr_abs.values, bins=timebinstot, histtype="step", label=lclabels[0])
for ite_sf in range(len(sf_times_init_goes)):
    ax_big.axvline(sf_times_init_goes[ite_sf], linestyle="--", color=colors[ite_sf], linewidth=1.0)
    ax_big.axvline(sf_times_end_goes[ite_sf], linestyle="--", color=colors[ite_sf], linewidth=1.0)
ax_big.set(xlabel="Time (s)", ylabel="count rate\n(couts/min)")
ax_big.legend()

# One solar flare in GOES per subplot
axes = []
for ite_sf in range(len(sf_times_init_goes)):
    row = 1 + int(ite_sf / 2)
    col = ite_sf % 2
    print(row, col)
    init_sftime, end_sf_time = sf_times_init_goes[ite_sf] - 10 * 60, sf_times_end_goes[ite_sf] + 10 * 60
    timebins = np.arange(init_sftime, end_sf_time+1, 1)
    ax = figfullsf.add_subplot(gs[row, col])
    ax.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=sf_names[ite_sf], color=colors[ite_sf])
    ax.axvline(sf_times_init_goes[ite_sf], linestyle="--", color=colors[ite_sf], linewidth=1.0)
    ax.axvline(sf_times_end_goes[ite_sf], linestyle="--", color=colors[ite_sf], linewidth=1.0)
    ax.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
    ax.legend()
    axes.append(ax)
plt.tight_layout(pad=0.8, h_pad=0.3, w_pad=0.8)
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/LC_solarflare_full.png")
plt.close(figfullsf)

# # GBM Solar flares
# sf_names = ["SFLARE24062467", "SFLARE24062458", "SFLARE24062426",
#             "SFLARE24062420", "SFLARE24062404", "SFLARE24062378",
#             "SFLARE24062348"]
# sf_dates = ["24/06/2024 16:09:14", "24/06/2024 13:56:15",
#             "24/06/2024 06:26:45", "24/06/2024 04:50:05",
#             "24/06/2024 01:05:37", "23/06/2024 18:45:53",
#             "23/06/2024 11:31:50"]
# # present in GOES ? [no, no, no, yes, no, no, yes]
# sf_times_abs = [166154, 158175, 131205, 125405, 111937, 89153, 63110]


# Plot of GOES_SF8 to define starting and ending times of the event in the data
figsf8, ax = plt.subplots(1, 1, figsize=(10, 6))
ite_sf = 7
init_sftime, end_sf_time = sf_times_init_goes[ite_sf] - 2 * 60, sf_times_end_goes[ite_sf] + 2 * 60
timebins = np.arange(init_sftime, end_sf_time+1, 1)
ax.hist(df_list[0].time_corr_abs.values, bins=timebins, histtype="step", label=sf_names[ite_sf], color=colors[ite_sf])
ax.axvline(sf_times_init_goes[ite_sf], linestyle="--", color=colors[ite_sf], linewidth=1.0)
ax.axvline(sf_times_end_goes[ite_sf], linestyle="--", color=colors[ite_sf], linewidth=1.0)
ax.axvline(150362, linestyle="--", color="red", linewidth=1.0)
ax.axvline(150416, linestyle="--", color="red", linewidth=1.0)
ax.axvline(150060, linestyle="--", color="green", linewidth=1.0)
ax.axvline(150150, linestyle="--", color="green", linewidth=1.0)
ax.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax.legend()
plt.tight_layout()
plt.show()

# df_sf8_bkg = df_list[0][(df_list[0].time_corr_abs >= 150060) & (df_list[0].time_corr_abs < 150139)]
# df_sf8_signal = df_list[0][(df_list[0].time_corr_abs >= 150362) & (df_list[0].time_corr_abs <= 150416)]
init_sfbkg, end_sfbkg = 150000, 150210
df_sf8_bkg = df_list[0][(df_list[0].time_corr_abs >= init_sfbkg) & (df_list[0].time_corr_abs < end_sfbkg)]
init_sfsignal, end_sfsignal = 150362, 150417
df_sf8_signal = df_list[0][(df_list[0].time_corr_abs >= init_sfsignal) & (df_list[0].time_corr_abs < end_sfsignal)]

numberbins = [30, 200]
for nbins in numberbins:
    # Spectrum background sf8
    energybins = np.linspace(15, 5500, nbins)
    figbkg_sf8, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.hist(df_sf8_bkg.energy.values, bins=energybins, histtype="step", label="Background spectrum", weights=[1 / len(df_sf8_bkg.energy.values)] * len(df_sf8_bkg.energy.values))
    ax.hist(df_sf8_signal.energy.values, bins=energybins, histtype="step", label="GOES_SF8 spectrum", weights=[1 / len(df_sf8_signal.energy.values)] * len(df_sf8_signal.energy.values))
    ax.axvline(511, linestyle="--", color="green", label="511 keV")
    ax.set(xlabel="Energy (keV)", ylabel="Normalised number of events", yscale="log")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/D2B_spec_sf8_nbins{nbins}.png")
    plt.close(figbkg_sf8)

    energybins = np.logspace(np.log10(15), np.log10(5500), nbins)
    figbkg_sf8_loglog, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.hist(df_sf8_bkg.energy.values, bins=energybins, histtype="step", label="Background spectrum", weights=[1 / len(df_sf8_bkg.energy.values)] * len(df_sf8_bkg.energy.values))
    ax.hist(df_sf8_signal.energy.values, bins=energybins, histtype="step", label="GOES_SF8 spectrum", weights=[1 / len(df_sf8_signal.energy.values)] * len(df_sf8_signal.energy.values))
    ax.axvline(511, linestyle="--", color="green", label="511 keV")
    ax.axvline(90, linestyle="--", color="orange", label="90 keV")
    ax.set(xlabel="Energy (keV)", ylabel="Normalised number of events", yscale="log", xscale="log")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/D2B_spec_sf8_nbins{nbins}_loglog.png")
    plt.close(figbkg_sf8_loglog)

    energybins = np.logspace(np.log10(15), np.log10(5500), nbins)
    figbkg_sf8_loglog_cr, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.hist(df_sf8_bkg.energy.values, bins=energybins, histtype="step", label="Background spectrum", weights=[1 / (end_sfbkg - init_sfbkg)] * len(df_sf8_bkg.energy.values))
    ax.hist(df_sf8_signal.energy.values, bins=energybins, histtype="step", label="GOES_SF8 spectrum", weights=[1 / (end_sfsignal - init_sfsignal)] * len(df_sf8_signal.energy.values))
    ax.axvline(511, linestyle="--", color="green", label="511 keV")
    # ax.axvline(90, linestyle="--", color="orange", label="90 keV")
    ax.set(xlabel="Energy (keV)", ylabel="Count rate (counts/s)", yscale="log", xscale="log")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/D2B_spec_sf8_nbins{nbins}_loglog_countrates.png")
    plt.close(figbkg_sf8_loglog_cr)

# df_sf8_bkg = df_list[0][(df_list[0].time_corr_abs >= 150060) & (df_list[0].time_corr_abs < 150139)]
# df_sf8_signal = df_list[0][(df_list[0].time_corr_abs >= 150362) & (df_list[0].time_corr_abs <= 150416)]
# print(len(df_sf8_signal) - len(df_sf8_bkg))
#
# energybins = np.linspace(1, 5500, 30)
# # Spectrum background sf8
# figbkg_sf8, ax = plt.subplots(1, 1, figsize=(10, 6))
# ax.hist(df_sf8_bkg.energy.values, bins=energybins, histtype="step", label="Background spectrum")
# ax.hist(df_sf8_signal.energy.values, bins=energybins, histtype="step", label="GOES_SF8 spectrum")
# ax.set(xlabel="Energy (keV)", ylabel="Number of events", yscale="log")
# ax.legend()
# plt.tight_layout()
# plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/spec/D2B_spec_sf8.png")
# plt.close(figbkg_sf8)



# TGF
tgfdata = "25/06/2024 23:53:42"
tgftime_abs = 280422
init_tgf_time, end_tgf_time = tgftime_abs - 1 * 60, tgftime_abs + 1 * 60
timebinstgf = np.arange(init_tgf_time, end_tgf_time + 1, 1)
plt.rcParams.update({'font.size': 15})
# Automatic all light curves
figlctgf, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True)
h1 = ax1.hist(df_list[0].time_corr_abs.values, bins=timebinstgf, histtype="step", label=lclabels[0])
ax1.set(ylabel="count rate\n(couts/s)")
ax1.legend()

ax2.hist(df_list[1].time_corr_abs.values, bins=timebinstgf, histtype="step", label=lclabels[1])
h2 = ax2.set(ylabel="count rate\n(couts/s)")
ax2.legend()

h3 = ax3.hist(np.concatenate((df_list[2].time_corr_abs.values, df_list[3].time_corr_abs.values, df_list[4].time_corr_abs.values, df_list[5].time_corr_abs.values)), bins=timebinstgf, histtype="step", label=f"D2A light curve")
ax3.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax3.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/triple_LC_tgf_{init_tgf_time}-{end_tgf_time}.png")
plt.close(figlctgf)

figlctgf, ax = plt.subplots(1, 1, figsize=(10, 3), sharex=True)
ax.hist(df_list[0].time_corr_abs.values, bins=timebinstgf, histtype="step", label=lclabels[0])
ax.set(xlabel="Time (s)", ylabel="count rate\n(couts/s)")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/LC/D2B_LC_tgf_{init_tgf_time}-{end_tgf_time}.png")
plt.close(figlctgf)



# Electron precipitation
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from datetime import datetime, timezone

mpl.use("Qt5Agg")
plt.rcParams.update({'font.size': 13})
filenames = ["corr_abs_maud_file.root", "corr_abs_dssd_file.root", "corr_abs_ucda_file.root",
             "corr_abs_ucdb_file.root", "corr_abs_ucdc_file.root", "corr_abs_ucdd_file.root"]
detec = ["D2B", "D1B", "D2AA", "D2AB", "D2AC", "D2AD"]
enerlabels = []
lclabels = []

df_list = []
begin = 3600
end = "end"
for fileite, filename in enumerate(filenames):
    with uproot.open(filename) as file:
        # here, file is a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        saved_columns = ["time_corr_abs", "pps_cpt_corr_abs", "energy"]

        # convert the tree in dataframe
        data = pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=100000, library="pd")], ignore_index=True)

    init_time = np.min(data["time_corr_abs"].values)
    final_time = np.max(data["time_corr_abs"].values)

    if begin == "begin":
        begin_time = int(init_time)
    else:
        begin_time = begin
    if end == "end":
        end_time = int(final_time) + 1
    else:
        end_time = end

    df_list.append(data[(data.time_corr_abs >= begin_time) & (data.time_corr_abs <= end_time)])
    enerlabels.append(f"{detec[fileite]} spectrum over {begin_time}-{end_time} time range")
    lclabels.append(f"{detec[fileite]} light curve")


ep_init, ep_end = 40000, 125000
ep_selec_init1, ep_selec_end1 = 57100, 70000
ep_selec_init2, ep_selec_end2 = 91000, 107000
ep_selec_init3, ep_selec_end3 = 110000, 122000
bkg1_init, bkg1_end = 42000, 55300
bkg2_init, bkg2_end = 77000, 85000
bkg3_init, bkg3_end = 86500, 90000
timebins = np.arange(ep_init, ep_end+60, 60)

ep_df = df_list[0][(df_list[0].time_corr_abs >= ep_init) & (df_list[0].time_corr_abs < ep_end)]
ep_selec_df1 = df_list[0][(df_list[0].time_corr_abs >= ep_selec_init1) & (df_list[0].time_corr_abs < ep_selec_end1)]
ep_selec_df2 = df_list[0][(df_list[0].time_corr_abs >= ep_selec_init2) & (df_list[0].time_corr_abs < ep_selec_end2)]
ep_selec_df3 = df_list[0][(df_list[0].time_corr_abs >= ep_selec_init3) & (df_list[0].time_corr_abs < ep_selec_end3)]
ep_bkg1_df = df_list[0][(df_list[0].time_corr_abs >= bkg1_init) & (df_list[0].time_corr_abs < bkg1_end)]
ep_bkg2_df = df_list[0][(df_list[0].time_corr_abs >= bkg2_init) & (df_list[0].time_corr_abs < bkg2_end)]
ep_bkg3_df = df_list[0][(df_list[0].time_corr_abs >= bkg3_init) & (df_list[0].time_corr_abs < bkg3_end)]

colors = ["blue", "orange", "green", "red", "purple", "pink", "gray", "olive", "cyan"]
figep, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(ep_df.time_corr_abs.values, bins=timebins, histtype="step", label=lclabels[0])
# EP 1
ax.axvspan(ep_selec_init1, ep_selec_end1, color=colors[0], alpha=0.3, label="EP1")
# EP 2
ax.axvspan(ep_selec_init2, ep_selec_end2, color=colors[1], alpha=0.3, label="EP2")
# EP 3
ax.axvspan(ep_selec_init3, ep_selec_end3, color=colors[2], alpha=0.3, label="EP3")
# BKG1
ax.axvspan(bkg1_init, bkg1_end, color=colors[3], alpha=0.3, label="Background1")
# BKG2
ax.axvspan(bkg2_init, bkg2_end, color=colors[4], alpha=0.3, label="Background2")
# BKG3
ax.axvspan(bkg3_init, bkg3_end, color=colors[5], alpha=0.3, label="Background3")
ax.set(xlabel="Time (s)", ylabel="count rate (couts/min)")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/elec_prep/LC_split.png")
plt.close(figep)
# plt.show()

energybins = np.linspace(15, 5500, 500)

figspecep, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(ep_selec_df1.energy.values, bins=energybins, histtype="step", label="EP1", weights=[1 / len(ep_selec_df1.energy.values)] * len(ep_selec_df1.energy.values), color=colors[0])
ax.hist(ep_selec_df2.energy.values, bins=energybins, histtype="step", label="EP2", weights=[1 / len(ep_selec_df2.energy.values)] * len(ep_selec_df2.energy.values), color=colors[1])
ax.hist(ep_selec_df3.energy.values, bins=energybins, histtype="step", label="EP3", weights=[1 / len(ep_selec_df3.energy.values)] * len(ep_selec_df3.energy.values), color=colors[2])
ax.hist(ep_bkg1_df.energy.values, bins=energybins, histtype="step", label="Background 1", weights=[1 / len(ep_bkg1_df.energy.values)] * len(ep_bkg1_df.energy.values), color=colors[3])
ax.hist(ep_bkg2_df.energy.values, bins=energybins, histtype="step", label="Background 2", weights=[1 / len(ep_bkg2_df.energy.values)] * len(ep_bkg2_df.energy.values), color=colors[4])
ax.hist(ep_bkg3_df.energy.values, bins=energybins, histtype="step", label="Background 3", weights=[1 / len(ep_bkg3_df.energy.values)] * len(ep_bkg3_df.energy.values), color=colors[5])
ax.set(xlabel="Energy (keV)", ylabel="Normalised number of events", yscale="log")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/elec_prep/D2B_ep_spec.png")
plt.close(figspecep)
# plt.show()

energybins = np.logspace(np.log10(15), np.log10(5500), 500)

figspeceplog, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(ep_selec_df1.energy.values, bins=energybins, histtype="step", label="EP1", weights=[1 / len(ep_selec_df1.energy.values)] * len(ep_selec_df1.energy.values), color=colors[0])
ax.hist(ep_selec_df2.energy.values, bins=energybins, histtype="step", label="EP2", weights=[1 / len(ep_selec_df2.energy.values)] * len(ep_selec_df2.energy.values), color=colors[1])
ax.hist(ep_selec_df3.energy.values, bins=energybins, histtype="step", label="EP3", weights=[1 / len(ep_selec_df3.energy.values)] * len(ep_selec_df3.energy.values), color=colors[2])
ax.hist(ep_bkg1_df.energy.values, bins=energybins, histtype="step", label="Background 1", weights=[1 / len(ep_bkg1_df.energy.values)] * len(ep_bkg1_df.energy.values), color=colors[3])
ax.hist(ep_bkg2_df.energy.values, bins=energybins, histtype="step", label="Background 2", weights=[1 / len(ep_bkg2_df.energy.values)] * len(ep_bkg2_df.energy.values), color=colors[4])
ax.hist(ep_bkg3_df.energy.values, bins=energybins, histtype="step", label="Background 3", weights=[1 / len(ep_bkg3_df.energy.values)] * len(ep_bkg3_df.energy.values), color=colors[5])
ax.axvline(511, linestyle="--", color="green", label="511 keV", linewidth=1.0)
ax.axvline(104, linestyle="--", color="red", label="104 keV", linewidth=1.0)
ax.set(xlabel="Energy (keV)", ylabel="Normalised number of events", xscale="log", yscale="log")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/elec_prep/D2B_ep_spec_loglog.png")
plt.close(figspeceplog)
# plt.show()

figspeceplogcr, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(ep_selec_df1.energy.values, bins=energybins, histtype="step", label="EP1", weights=[1 / (ep_selec_end1 - ep_selec_init1)] * len(ep_selec_df1.energy.values), color=colors[0])
ax.hist(ep_selec_df2.energy.values, bins=energybins, histtype="step", label="EP2", weights=[1 / (ep_selec_end2 - ep_selec_init2)] * len(ep_selec_df2.energy.values), color=colors[1])
ax.hist(ep_selec_df3.energy.values, bins=energybins, histtype="step", label="EP3", weights=[1 / (ep_selec_end3 - ep_selec_init3)] * len(ep_selec_df3.energy.values), color=colors[2])
ax.hist(ep_bkg1_df.energy.values, bins=energybins, histtype="step", label="Background 1", weights=[1 / (bkg1_end - bkg1_init)] * len(ep_bkg1_df.energy.values), color=colors[3])
ax.hist(ep_bkg2_df.energy.values, bins=energybins, histtype="step", label="Background 2", weights=[1 / (bkg2_end - bkg2_init)] * len(ep_bkg2_df.energy.values), color=colors[4])
ax.hist(ep_bkg3_df.energy.values, bins=energybins, histtype="step", label="Background 3", weights=[1 / (bkg3_end - bkg3_init)] * len(ep_bkg3_df.energy.values), color=colors[5])
ax.axvline(511, linestyle="--", color="green", label="511 keV", linewidth=1.0)
# ax.axvline(104, linestyle="--", color="red", label="104 keV", linewidth=1.0)
ax.set(xlabel="Energy (keV)", ylabel="Count rate (counts/s)", xscale="log", yscale="log")
ax.legend()
plt.tight_layout()
plt.savefig(f"/home/nathan/Desktop/GRB_polarimetry/Thesis/reports_kiruna/results/elec_prep/D2B_ep_spec_loglog_countrates.png")
plt.close(figspeceplogcr)
# plt.show()
