import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import rfft, rfftfreq
import uproot
import pandas as pd
from scipy.stats import poisson
import matplotlib as mpl
from memory_profiler import memory_usage
from datetime import datetime
from scipy.optimize import curve_fit

mpl.use('TkAgg')


############################################################################################################
# For information on uproot :
# https://uproot.readthedocs.io/en/latest/basic.html
############################################################################################################

class KirunaAnalysis:
    """
    Class containing Kiruna data and method to easily manipulate them
    """
    def __init__(self, d1a_datafile=None, d2a_datafile=None, d1b_datafile=None, d2b_datafile=None):
        """

        :param d2b_datafile:
        """
        self.epoch_ref = 1719079200
        # Origin of the time referential : 2024/06/22 20:00:00
        self.timeref = datetime(2024, 6, 22, 20, 0, 0).timestamp()

        # CEA data
        self.d1a_energy_cutted = None
        if d1a_datafile is None:
            self.d1a_datafile = ""
        else:
            self.d1a_datafile = d1a_datafile
        self.d1a_df = None
        self.d1a_init_ts = None
        self.d1a_final_ts = None
        self.d1a_init_date = None
        self.d1a_final_date = None

        # DSSD data
        self.d1b_energy_cutted = None
        if d1b_datafile is None:
            self.d1b_datafile = "corr_abs_dssd_file.root"
        else:
            self.d1b_datafile = d1b_datafile
        self.d1b_df = None
        self.d1b_init_ts = None
        self.d1b_final_ts = None
        self.d1b_init_date = None
        self.d1b_final_date = None

        # UCD data
        self.d2a_energy_cutted = None
        if d2a_datafile is None:
            self.d2a_datafile = ["corr_abs_ucda_file.root", "corr_abs_ucdb_file.root", "corr_abs_ucdc_file.root", "corr_abs_ucdd_file.root"]
        else:
            self.d2a_datafile = d2a_datafile
        self.d2a_df = None
        self.d2a_init_ts = None
        self.d2a_final_ts = None
        self.d2a_init_date = None
        self.d2a_final_date = None

        # Maud data
        self.d2b_energy_cutted = None
        if d2b_datafile is None:
            self.d2b_datafile = "corr_abs_maud_file.root"
        else:
            self.d2b_datafile = d2b_datafile
        self.d2b_df = None
        self.d2b_init_ts = None
        self.d2b_final_ts = None
        self.d2b_init_date = None
        self.d2b_final_date = None

        # Mean flux of 0.3 ph/s but emission every 0.03 sec : mean num of ph emited per pulsation : 0.3*0.03 = 0.009
        self.crab_flux = 0.02
        self.crab_period = 33.9e-3
        self.crab_frequency = round(1 / self.crab_period, 3)
        # Total flux between 320000-323000 s : 506449
        # Source flux of 0.3 ph/s so 900 ph over this duration
        # >>> background flux integrated over 3000s : 505549  >>>> mean flux of 168.52 ph/s but the bin is not 1s wide so timescale corrects
        self.bkg_flux_3000s = 168.52

        self.fft_binning = 1e-2
        self.fft_signal = None
        self.emptyness = [[19237.272103425, 21130.70923795], [22409.7524636, 22434.4893218], [26444.113075475, 28864.92999595],
                         [64440.98104165, 64441.6488771], [85554.701873725, 85891.042379775], [248483.597786375, 248630.742785725],
                         [248800.05419025, 248883.775494775], [249157.926650425, 249382.772406325], [251492.077928625, 251649.213788825],
                         [260324.563231225, 260327.554940825], [261335.77597295, 261336.297486475], [323258.856362225, 323304.646581125]]

        self.make_d2a_df()
        # self.make_d2b_df()

    def make_d2a_df(self):
        """
        Extract the d2b panda df from a rootfile
        :return:
        """
        # saved_columns = ["pps_cpt_corr_abs", "time_corr_abs", "energy", "ts_init", "pps_cpt_init", "adc_channels"]
        saved_columns = ["pps_cpt_corr_abs", "time_corr_abs", "energy", "pps_cpt_init"]
        temp_df = []
        for udcfile in self.d2a_datafile:
            with uproot.open(udcfile) as file:
                print("Memory 1 : ", memory_usage(), " MB")
                tree = file["Events"]
                print("Memory 2 : ", memory_usage(), " MB")

                # convert the tree in dataframe
                temp_df.append(pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=1000000, library="pd")], ignore_index=True))
            print("Memory 3 : ", memory_usage(), " MB")
        subdets = ["ucda", "ucdb", "ucdc", "ucdd"]
        print(len(temp_df))
        for df_ite, df in enumerate(temp_df):
            df["sub_det"] = subdets[df_ite]
            print(f"subdet : {subdets[df_ite]} - len : {len(df)}")
        # self.d2a_df = pd.concat(temp_df, ignore_index=True).sort_values(by=["time_corr_abs", "sub_det"], inplace=True)
        self.d2a_df = pd.concat(temp_df, ignore_index=True)
        print(self.d2a_df)
        self.d2a_df.sort_values(by=["time_corr_abs", "sub_det"], inplace=True)
        print(self.d2a_df)
        self.d2a_df.reset_index(drop=True, inplace=True)
        print(self.d2a_df)
        print(len(self.d2a_df))
        self.d2a_init_ts = self.d2a_df.pps_cpt_corr_abs.values[0]
        self.d2a_final_ts = self.d2a_df.pps_cpt_corr_abs.values[-1]
        self.d2a_init_date = datetime.fromtimestamp(self.d2a_init_ts + self.timeref)
        self.d2a_final_date = datetime.fromtimestamp(self.d2a_final_ts + self.timeref)

    def make_d2b_df(self):
        """
        Extract the d2b panda df from a rootfile
        :return:
        """
        # saved_columns = ["pps_cpt_corr_abs", "time_corr_abs", "energy", "ts_init", "pps_cpt_init", "adc_channels"]
        saved_columns = ["pps_cpt_corr_abs", "time_corr_abs", "energy", "ts_init", "pps_cpt_init", "energy_adc"]
        with uproot.open(self.d2b_datafile) as file:
            # file is there a TDirectory, you can display the content with :
            # print(file.keys())
            # You can have even more detail on what the key is with :
            # print(file.classnames())
            # Accessing the tree
            print("Memory 1 : ", memory_usage(), " MB")
            tree = file["Events"]
            print("Memory 2 : ", memory_usage(), " MB")
            # Maud tree keys :
            # print(maud_tree.keys())

            # convert the tree in dataframe
            self.d2b_df = pd.concat([chunk for chunk in tree.iterate(saved_columns, step_size=1000000, library="pd")], ignore_index=True)
            # self.d2b_df = tree.arrays(library="pd", entry_start=1609950, entry_stop=1609955)
            print("Memory 3 : ", memory_usage(), " MB")

            self.d2b_init_ts = self.d2b_df.pps_cpt_corr_abs.values[0]
            self.d2b_final_ts = self.d2b_df.pps_cpt_corr_abs.values[-1]
            self.d2b_init_date = datetime.fromtimestamp(self.d2b_init_ts + self.timeref)
            self.d2b_final_date = datetime.fromtimestamp(self.d2b_final_ts + self.timeref)

    def energy_cut(self, erg_cut):
        """
        Performs an energy cut on the data
        :param erg_cut:
        :return:
        """
        if self.energy_cutted is not None:
            if erg_cut[0] < self.energy_cutted[0] or erg_cut[1] > self.energy_cutted[1]:
                raise ValueError("Error during the energy cut : a cut has already been applied, the new cut is wider than the previous one so data is incomplete")
            elif erg_cut[0] == self.energy_cutted[0] or erg_cut[1] == self.energy_cutted[1]:
                print("This energy cut is already applied, values are kept as they were")
                return
            else:
                print("Warning : a cut has already been applied but the 2 cuts should not interfere, the data should be complete")
        self.energy_cutted = erg_cut
        self.d2b_df = self.d2b_df[np.logical_and(self.d2b_df.energy > erg_cut[0], self.d2b_df.energy < erg_cut[1])]

    def fft_analysis(self, erg_cut=None, mode="kiruna", add_signal=False, times_of_emptyness=False, periods=None, figtitle=""):
        """

        :param erg_cut:
        :param mode:
        :param add_signal:
        :param times_of_emptyness:
        :param periods:
        :param figtitle:
        :return:
        """
        if erg_cut is not None:
            self.energy_cut(erg_cut)
            figtitle = f"{figtitle} - ergcut : {erg_cut}"
        else:
            figtitle = f"{figtitle} - ergcut : None"
        if periods is None:
            periods = [["23/06/2024 4:20:00", "23/06/2024 18:10:00"],
                       ["24/06/2024 5:00:00", "24/06/2024 21:00:00"],
                       ["25/06/2024 6:40:00", "25/06/2024 22:50:00"],
                       ["26/06/2024 8:00:00", "end"]]
        deg_lim = 15
        kept_times = [[self.get_timestamp(datestring) for datestring in plist] for plist in periods]

        if mode == "kiruna":
            # Opening ROOT file
            name = f"{figtitle}\nFFT with corrected time"
            self.get_bined_times(name, add_signal=add_signal, emptyness=times_of_emptyness)
            self.run_fft(name, vline=True)
        elif mode == "simulated":
            name = f"{figtitle}\nFFT with fake signal"
            self.get_fake_signal(name, add_signal=add_signal, emptyness=times_of_emptyness)
            self.run_fft(name, vline=True)
        elif mode == "kiruna_select":
            name = f"{figtitle}\nFFT with corrected time"
            self.get_bined_times(name, add_signal=add_signal, emptyness=times_of_emptyness, time_ranges=kept_times)
            self.run_fft(name, namefig=f"window_{len(periods)}_{deg_lim}°", vline=True)
        else:
            raise ValueError("Please use a correct mode for the fft analysis : kiruna, kiruna_select or simulated")

    def create_spectrum(self, date_init, date_final, erg_cut=None, nbins=None, xlog=False, ylog=False, fit511=True):
        """

        :param date_init:
        :param date_final:
        :param erg_cut:
        :param xlog:
        :param ylog:
        :param fit511:
        :return:
        """
        if erg_cut is not None:
            self.energy_cut(erg_cut)
        init_date = self.get_timestamp(date_init)
        finish_date = self.get_timestamp(date_final)

        if xlog:
            xlog = "log"
        else:
            xlog = "linear"
        if ylog:
            ylog = "log"
        else:
            ylog = "linear"

        maud_df_select = self.d2b_df[np.logical_and(self.d2b_df.pps_cpt_corr_abs >= init_date, self.d2b_df.pps_cpt_corr_abs <= finish_date)]
        if len(maud_df_select) > 0:
            # ============================================================================================================================
            # Ploting the spectrum
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))

            if fit511:
                if nbins is not None:
                    bins = np.linspace(-100, 5500, nbins)
                else:
                    bins = np.linspace(-100, 5500, 1000)
                yhist, xhist = np.histogram(maud_df_select.energy.values, bins=bins)

                ener_select, popt, r2 = fitting511(xhist, yhist, 511, 430, 590)
                print(f"============ Fit with {len(bins)} bins ============")
                print("===  Fitting results")
                if len(popt) == 0:
                    print(f"   Mean = {0} keV       Energy resolution = {0} %  ==")
                else:
                    print(f"   Mean = {popt[1]} keV       Energy resolution = {popt[2] * 2.3548 / popt[1] * 100} %  ==")
                    ax.scatter(ener_select, gauss_func(ener_select, *popt), color="red", s=2)
                print(f"   R² = {r2}")
                # print("===  Positions")
                # print(f"   init = {valinit} keV")
                # print(f"   511 = {val511} keV")
                # print(f"   end = {valfinal} keV")
                print()
                ax.axvline(511)
                ax.axvline(430)
                ax.axvline(590)
                # for val in change:
                #     ax.axvline(val)
            else:
                bins = np.linspace(-100, 5500, 1000)
            ax.hist(maud_df_select.energy.values, bins=bins, histtype="step", label="In flight spectrum")
            ax.axvline(511, label="511 keV", color="green")
            ax.set(xlabel="Energy (keV)", ylabel="Number of event", xscale=xlog, yscale=ylog)
            ax.legend()
            plt.show()
        else:
            print("No data for the given period")

    def create_lightcurve(self, init_date="begin", end_date="end", erg_cut=None, bins=1000):
        """

        :param init_date:
        :param end_date:
        :param erg_cut:
        :return:
        """
        if erg_cut is not None:
            self.energy_cut(erg_cut)

        init_date = self.get_timestamp(init_date)
        finish_date = self.get_timestamp(end_date)

        maud_df_select = self.d2b_df[np.logical_and(self.d2b_df.pps_cpt_corr_abs >= init_date, self.d2b_df.pps_cpt_corr_abs <= finish_date)]
        if len(maud_df_select) > 0:
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            ax.hist(maud_df_select.time_corr_abs.values, bins=bins, histtype="step", label=f"In flight light curve\nergcut = {erg_cut}")
            # ax.axvline(511, label="511 keV", color="green")
            ax.set(xlabel="Corrected time from 2024/06/22 20:00:00 (s)", ylabel="Number of event")
            # ax.scatter(ener2_select, gauss_func(ener2_select, *popt2), color="orange", s=2)

            ax.legend()
            plt.show()

    def fit_511(self, graphs=True):
        """

        :return:
        """
        showfit = False

        r2lim = 0.75
        errorlim = 5
        nev_min = 30000
        timestep = 5 * 60

        saved_ratio = []
        saved_date = []
        saved_r2 = []
        saved_error = []
        saved_nev = []
        saved_reso = []

        loopcount = 0
        loopunfit = 0
        oldlen = 0
        init_ts = self.d2b_init_ts
        end_ts = self.d2b_init_ts + timestep
        # init_ts = self.d2b_init_ts + 90*timestep
        # end_ts = self.d2b_init_ts + 91*timestep

        # finish_ts = init_ts + 10*timestep
        finish_ts = self.d2b_final_ts
        while end_ts <= finish_ts and loopcount < 2000:
            maud_df_select = self.d2b_df[np.logical_and(self.d2b_df.pps_cpt_corr_abs >= init_ts, self.d2b_df.pps_cpt_corr_abs <= end_ts)]
            nevents = len(maud_df_select)
            lr2, lerror = [], []
            if nevents <= nev_min:
                ener_select, popt, r2, bins, val511, valinit, valfinal, change = [], [], 0, [], 0, 0, 0, []
                print(f"Not enough data between {init_ts} and {end_ts} : no fit done")
            elif loopunfit != 0 and nevents == oldlen:
                ener_select, popt, r2, bins, val511, valinit, valfinal, change = [], [], 0, [], 0, 0, 0, []
                # print("Making the sample bigger didn't change the number of values - fit is not performed (added ts has no events)")
            else:
                print(f"=== {nevents} events for the time period {init_ts} - {end_ts}")
                lener_select, lpopt, lr2, lerror, lbins, l511, linit, lfinal, lchange = findbestfit(maud_df_select.energy.values, detailed=True)
                idmaxr2 = np.argmax(np.where(np.array(lerror) < errorlim, np.array(lr2), 0))
                # print(idmaxr2)
                if showfit and fit_validation(lr2, lerror, r2lim, errorlim):
                    # ============================================================================================================================
                    # Ploting the spectrum
                    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 10))
                    fig.suptitle(f"Period {init_ts} - {end_ts}")
                    xlog = "linear"
                    ylog = "linear"
                    axes = [ax1, ax2, ax3, ax4]
                    tested_nbins = [900, 700, 500, 300]
                    for ite in range(len(lbins)):
                        print(f"============ Fit with {len(lbins[ite])} bins ============")
                        print("===  Fitting results")
                        if len(lpopt[ite]) == 0:
                            print(f"   Mean = {0} keV       Energy resolution = {0} %  ==")
                        else:
                            print(f"   Mean = {lpopt[ite][1]} keV       Energy resolution = {lpopt[ite][2] * 2.3548 / lpopt[ite][1] * 100} %  ==")
                            print(f"       Expected - Real ratio : {round(511 / lpopt[ite][1], 3)}")
                        print(f"   R² = {lr2[ite]}")
                        print("===  Positions")
                        print(f"   init = {linit[ite]} keV")
                        print(f"   511 = {l511[ite]} keV")
                        print(f"   end = {lfinal[ite]} keV")
                        print()

                        if len(lpopt[ite]) > 0:
                            axes[ite].axvline(l511[ite], color="violet", label=f"511 maximum : {round(l511[ite], 3)} keV")
                            axes[ite].axvline(lpopt[ite][1], color="cyan", label=f"511 fitted : {round(lpopt[ite][1], 3)} keV")
                            axes[ite].axvline(linit[ite], color="green", label="511 fit interval edges")
                            axes[ite].axvline(lfinal[ite], color="green")
                            # for val in lchange[ite]:
                            #     axes[ite].axvline(val, color="orange", linestyle='--')
                            # axes[ite].axvline(lchange[ite][0], color="orange", linestyle='--', label="minimum")
                            # axes[ite].axvline(lchange[ite][1], color="red", linestyle='--', label="maximum")
                            axes[ite].scatter(lener_select[ite], gauss_func(lener_select[ite], *lpopt[ite]), color="red", s=2, label=f"Fit : R² = {round(lr2[ite], 3)}\n      511peak-511 relative error = {round(lerror[ite], 3)} %")
                            if fit_validation([lr2[ite]], [lerror[ite]], r2lim, errorlim):
                                axes[ite].hist(maud_df_select.energy.values, bins=lbins[ite], histtype="step", color="green", label="In flight spectrum")
                            else:
                                axes[ite].hist(maud_df_select.energy.values, bins=lbins[ite], histtype="step", color="blue", label="In flight spectrum")
                            axes[ite].axvline(511, label="511 keV", color="black")
                            axes[ite].legend()
                        axes[ite].set(xlabel="Energy (keV)", ylabel="Number of event", xscale=xlog, yscale=ylog, title=f"Fit with nbins = {tested_nbins[ite]}")
                    axes[idmaxr2].text(0.5, 0.9, 'Selected fit', fontsize=12, color='red', alpha=0.7, ha='center', va='bottom', transform=axes[idmaxr2].transAxes)
                    plt.show()
            if fit_validation(lr2, lerror, r2lim, errorlim):  # Possibly another info to see if the fit is really ok !
                ratio = 511 / lpopt[idmaxr2][1]
                if len(saved_ratio) > 0:
                    old_ratio = saved_ratio[-1]
                else:
                    old_ratio = ratio
                if np.abs((old_ratio - ratio) / old_ratio) < 0.1:
                    saved_ratio.append(ratio)
                    saved_date.append(f"{init_ts} and {end_ts}")
                    saved_r2.append(lr2[idmaxr2])
                    saved_error.append(lerror[idmaxr2])
                    saved_nev.append(nevents)
                    saved_reso.append(lpopt[idmaxr2][2] * 2.3548 / lpopt[idmaxr2][1] * 100)
                    init_ts = end_ts
                    end_ts += timestep
                    loopcount += 1
                    loopunfit = 0
                    oldlen = 0
                    print(f"Saving value between {init_ts} and {end_ts} ")
                else:
                    if loopunfit == 0:
                        print(f"No fit with proper ratio for value between {init_ts} and {end_ts}  -  fine iteration to find a correct fit")
                    end_ts += 1
                    loopunfit += 1
                    oldlen = nevents

            else:
                if loopunfit == 0:
                    print(f"No fit for value between {init_ts} and {end_ts}  -  fine iteration to find a correct fit")
                end_ts += 1
                loopunfit += 1
                oldlen = nevents

        if end_ts != finish_ts:
            print(f"Saving value of this range as the one of previous period ?")

        if graphs:
            fit_time = (np.array([int(date.split(" ")[-1]) for date in saved_date]) + np.array([int(date.split(" ")[0]) for date in saved_date])) / 2

            fig1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex=True)
            ax1.axhline(1, color="red")
            ax1.plot(fit_time, saved_ratio, color="blue", label="Original")
            # ax1.scatter(fit_time, save_ratio, color="pink", s=10)
            # ax1.plot(fit_time, no_outliers, color="green")
            # ax1.plot(fit_time, smoothed3, color="orange")
            # ax1.scatter(fit_time, smoothed3, color="pink", s=10)
            # ax1.plot(fit_time, no_outliers_smoothed3, color="red", label="Smoothed 3 points")
            # ax1.plot(fit_time, no_outliers_smoothed5, color="orange", label="Smoothed 5 points")
            # ax1.scatter(range(init_ts, final_ts + 1), interpol_ratio, color="black", s=5)
            ax1.legend()
            ax1.set(ylabel="511 ratio")
            ax2.plot(fit_time, saved_reso, color="blue")
            ax2.set(ylabel="Fit resolution (%)")
            ax3.plot(fit_time, saved_nev)
            ax3.set(ylabel="N events")
            ax4.plot(fit_time, saved_r2)
            ax4.set(ylabel="R²")
            ax5.plot(fit_time, saved_error)
            ax5.set(xlabel="Time (s)", ylabel="Relative error (%)")
            plt.show()

            fig2, (ax11, ax22, ax33) = plt.subplots(3, 1)
            ax11.axhline(1, color="red")
            ax11.plot(fit_time, saved_ratio, color="blue", label="Original")
            ax11.legend()
            ax11.set(ylabel="511 ratio")
            ax22.plot(fit_time, saved_reso, color="blue")
            ax22.set(xlabel="Time (s)", ylabel="Fit resolution (%)")

            ax33.hist(saved_reso, bins=30)
            ax33.set(xlabel="Fit resolution (%)", ylabel="Count")

            plt.show()

        print("MEAN RESOLUTION : ", np.mean(saved_reso))
        return saved_ratio, saved_date, saved_r2, saved_error, saved_nev, saved_reso


#                 Enregistre qu il n y a pas de correction ici
#
# Si le nombre d'ev ne change pas apres avoir itere parce qu'on avait pas de fit alors il faut prendre la valeur trouvée juste avant
# puis reprendre a la premiere seconde qui a de nouveau des valeurs !
    # solutions :
    # Augmenter temps
    # faire des fits avec differents bins
    # On a peut etre que besoin du 511 et mettre un nombre de sigma autour ? sinon utiliser le minimum trouvé mais il a moins de sens avec moins de valeurs
    # Ameliorer recherche de minimum ?....

    def get_bined_times(self, fft_title, add_signal=False, emptyness=False, time_ranges=None):
        """

        :param fft_title:
        :param add_signal:
        :param emptyness:
        :param time_ranges:
        :return:
        """
        time_array = self.d2b_df.time_corr_abs.values

        print("===============================================================================================================")
        print(f" Preparing the fft : {fft_title}")
        init_time = round(int(time_array[0] / self.fft_binning) * self.fft_binning, 3)
        end_time = round((int(time_array[-1] / self.fft_binning) + 1) * self.fft_binning, 3)
        n_val = int((end_time - init_time) / self.fft_binning + 1)
        bins = np.linspace(init_time, end_time, n_val)
        print(f"  == init time      : {init_time}")
        print(f"  == end_time       : {end_time}")
        print(f"  == number of bins : {n_val - 1}")

        if time_ranges is not None:
            visible_array = []
            for val in time_array:
                seen = False
                for time_range in time_ranges:
                    if time_range[0] <= val <= time_range[1]:
                        seen = True
                if seen:
                    visible_array.append(val)
            visible_array = np.array(visible_array)
            data_signal = np.histogram(visible_array, bins=bins)[0]
        else:
            data_signal = np.histogram(time_array, bins=bins)[0]

        if add_signal:
            centroids = np.around(bins[:-1] % self.crab_period, 5)
            added_signal = self.crab_draw(centroids, self.crab_flux, self.fft_binning)
            if emptyness:
                correction_list = []
                for ite in range(len(bins) - 1):
                    verif_bool = False
                    for inf_lim, sup_lim in self.emptyness:
                        if bins[ite] >= inf_lim and bins[ite + 1] <= sup_lim:
                            verif_bool = True
                    if verif_bool:
                        correction_list.append(0)
                    else:
                        correction_list.append(1)
                added_signal = added_signal * np.array(correction_list)
            self.fft_signal = data_signal + added_signal
        else:
            self.fft_signal = data_signal

    def get_fake_signal(self, fft_title, add_signal=False, emptyness=False):
        """

        :param fft_title:
        :param add_signal:
        :param emptyness:
        :return:
        """
        print("===============================================================================================================")
        print(f" Preparing the fft : {fft_title}")
        init_time = self.d2b_df.time_corr_abs.values[0]
        end_time = self.d2b_df.time_corr_abs.values[-1]
        n_val = int((end_time - init_time) / self.fft_binning + 1)
        bins = np.linspace(init_time, end_time, n_val)
        print(f"  == init time      : {init_time}")
        print(f"  == end_time       : {end_time}")
        print(f"  == number of bins : {n_val - 1}")

        centroids = np.around(bins[:-1] % self.crab_period, 5)
        bkgflux = self.bkg_flux_3000s * self.fft_binning
        bkg_signal = np.array([poisson.rvs(bkgflux) for ite in range(len(bins) - 1)])
        if add_signal:
            # Mean flux of 0.3 ph/s but emission every 0.03 sec : mean num of ph emited per pulsation : 0.3*0.03 = 0.009
            added_signal = self.crab_draw(centroids, self.crab_flux, self.fft_binning)
            bkg_signal = bkg_signal + added_signal
        if emptyness:
            correction_list = []
            for ite in range(len(bins) - 1):
                verif_bool = False
                for inf_lim, sup_lim in self.emptyness:
                    if bins[ite] >= inf_lim and bins[ite + 1] <= sup_lim:
                        verif_bool = True
                if verif_bool:
                    correction_list.append(0)
                else:
                    correction_list.append(1)
            self.fft_signal = bkg_signal * np.array(correction_list)
        else:
            self.fft_signal = bkg_signal

    def crab_draw(self, xlist, mean, bin_timescale):
        """

        :param xlist:
        :param mean:
        :param bin_timescale:
        :return:
        """
        added_signal = []
        for centroid in xlist:
            if centroid < bin_timescale/2:
                added_signal.append(poisson.rvs(mean))
            elif centroid > self.crab_period-bin_timescale/2:
                added_signal.append(poisson.rvs(mean))
            else:
                added_signal.append(0)
        added_signal = np.array(added_signal)
        return added_signal

    def run_fft(self, fft_title, namefig=None, vline=True):
        """

        :param fft_title:
        :param namefig:
        :param vline:
        :return:
        """
        if self.fft_signal is None:
            raise ValueError("There is no signal to analyse")
        print(f" Running a Fourier Fast Transform analysis : {fft_title}")
        # Calcul de la TF rapide (FFT)
        func_fftvals = rfft(self.fft_signal)

        x_axis = rfftfreq(len(self.fft_signal), self.fft_binning)
        y_axis = np.abs(func_fftvals)
        index = np.where(x_axis > 0.1, True, False)

        new_x, new_y = x_axis[index], y_axis[index]

        title = f"{fft_title}, bining timescale : {self.fft_binning}"
        # Visualisation de la TF
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        if vline:
            ax.axvline(self.crab_frequency, color="red", label="Crab pulsar frequency")
        ax.plot(new_x, new_y, label="FFT of the signal")
        ax.set(xlabel='Frequency [Hz]', ylabel='Amplitude', title=title)
        ax.legend()
        if namefig is not None:
            plt.savefig(namefig)
        plt.show()

    def get_timestamp(self, date):
        """

        :param date:
        :return:
        """
        if date == "begin":
            return self.d2b_df.pps_cpt_corr_abs.values[0]
        elif date == "end":
            return self.d2b_df.pps_cpt_corr_abs.values[-1]
        else:
            dateval, timeval = date.split(" ")
            day, month, year = map(int, dateval.split("/"))
            hour, mins, sec = map(int, timeval.split(":"))
            return datetime(year, month, day, hour, mins, sec).timestamp() - self.timeref

    def check_data_gaps(self):
        """

        :return:
        """
        time_array = self.d2b_df.time_corr_abs.values
        ite_list = []
        for ite in range(len(time_array) - 1):
            if time_array[ite + 1] - time_array[ite] > 0.5:
                ite_list.append([time_array[ite], time_array[ite + 1]])
        print(ite_list)
        check_times = False
        if check_times:
            gps_from_pps_array = self.d2b_df.pps_cpt_corr_abs.values + self.epoch_ref
            diff_list = []
            for ite in range(len(time_array) - 1):
                if time_array[ite + 1] - time_array[ite] > 1:
                    print(f"Diff : {time_array[ite + 1] - time_array[ite]}   ite {ite}-{ite + 1}/{len(time_array)}")
                    diff_list.append(ite)
                    date_time_local = datetime.fromtimestamp(gps_from_pps_array[ite])
                    date_time_localp1 = datetime.fromtimestamp(gps_from_pps_array[ite + 1])
                    print(f"Date et heure locale : {date_time_local} - {date_time_localp1}\n")
                if time_array[ite + 1] <= time_array[ite]:
                    print("ERROR")

            print(f"Num of vals in a row : {diff_list[0]}")
            for ite in range(len(diff_list) - 1):
                print(f"Num of vals in a row : {diff_list[ite + 1] - diff_list[ite]}")
            print(f"Num of vals in a row : {len(time_array) - diff_list[-1]}")


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
    return ampl * np.exp(-0.5*((x-mu)/sigma)**2) + a * (x - x0) + b


def fitting511(xhist, yhist, val511, lower_x, higher_x):
    """

    :param xhist:
    :param yhist:
    :param lower_x:
    :param higher_x:
    :return:
    """
    idxmin = np.where(xhist > lower_x)[0][0]
    adxmax = np.where(xhist < higher_x)[0][-1]

    x_select, y_select = xhist[max(idxmin - 1, 0):min(adxmax + 1, len(xhist))], yhist[max(idxmin - 1, 0):adxmax]
    x_select = (x_select[1:] + x_select[:-1]) / 2
    # print("number : ", np.sum(y_select))
    # print(val511, lower_x, higher_x)
    # print(idxmin, adxmax, f"{idxmin - 1}:{adxmax + 1}", f"{idxmin - 1}:{adxmax}")
    # print()
    a1, b1, x01 = (y_select[-1] - y_select[0]) / (x_select[-1] - x_select[0]), y_select[0], x_select[0]
    if a1 >= 0:
        return [], [], 0
    else:
        inf_bounds = [0, lower_x, 0, a1 * 1.1, b1 * 0.9, x01 * 0.9]
        sup_bounds = [np.inf, higher_x, val511-lower_x, a1 * 0.9, b1 * 1.1, x01 * 1.1]
        # guess = [1000, val511, 15, a1, b1, x01]
        guess = [10 * np.sum(y_select) / 1500, val511, 15, a1, b1, x01]

        try:
            popt, pcov = curve_fit(gauss_func, x_select, y_select, bounds=(inf_bounds, sup_bounds), p0=guess)[:2]
            y_fitted = gauss_func(x_select, *popt)

            # residual sum of squares
            ss_res = np.sum((y_select - y_fitted) ** 2)
            # total sum of squares
            ss_tot = np.sum((y_select - np.mean(y_select)) ** 2)
            # r-squared
            r2 = 1 - (ss_res / ss_tot)
            # print(f"R² = {r2}")
            # for itepar in range(len(popt)):
            #     print(f"{itepar} : {popt[itepar]} +- {np.sqrt(pcov[itepar][itepar])}")
        except RuntimeError:
            print("The fit does not converge")
            popt, pcov, r2 = [], [], 0

        return x_select, popt, r2


def findbestfit(energies, detailed=False):
    """

    :param energies:
    :param detailed:
    :return:
    """
    lbins = [np.linspace(-100, 5500, nbin) for nbin in [1000, 750, 500, 250]]
    lchange = []
    lchangeite = []
    linit = []
    l511 = []
    lfinal = []
    lener_select = []
    lpopt = []
    lr2 = []
    lerror = []
    for bins in lbins:
        yhist, xhist = np.histogram(energies, bins=bins)
        countp, enerp = yhist[1:] - yhist[:-1], (xhist[1:] + xhist[:-1]) / 2
        change, changeite = local_optima_finder(countp, enerp, yhist, "adapted")

        lchange.append(change)
        lchangeite.append(changeite)
        if len(change) == 2:
            peak511, minmaxdist = enerp[changeite[1]], enerp[changeite[1]] - enerp[changeite[0]]
            init_ener = peak511 - 1.4 * minmaxdist
            final_ener = peak511 + 1.4 * minmaxdist
            ener_select, popt, r2 = fitting511(xhist, yhist, peak511, init_ener, final_ener)
            linit.append(init_ener)
            l511.append(peak511)
            lfinal.append(final_ener)
            if r2 > 0:
                lener_select.append(ener_select)
                lpopt.append(popt)
                lr2.append(r2)
                lerror.append(np.abs((peak511 - popt[1]) / peak511) * 100)
            else:
                lener_select.append([])
                lpopt.append([])
                lr2.append(0)
                lerror.append(100)
        else:
            linit.append(0)
            l511.append(0)
            lfinal.append(0)
            lener_select.append([])
            lpopt.append([])
            lr2.append(0)
            lerror.append(100)
    if detailed:
        return lener_select, lpopt, lr2, lerror, lbins, l511, linit, lfinal, lchange
    else:
        best_index = np.argmax(lr2)
        return lener_select[best_index], lpopt[best_index], lr2[best_index], lerror[best_index], lbins[best_index], l511[best_index], linit[best_index], lfinal[best_index], lchange[best_index]


def local_optima_finder(der, ener, count, mode):
    """

    :param der:
    :param ener:
    :param count:
    :param mode:
    :return:
    """
    change = []
    changeite = []
    if mode == "global":
        for ite in range(1, len(der[:-1])):
            if ite < len(der) - 3:
                if der[ite] * der[ite + 1] < 0:
                    if 300 < ener[ite] < 1000:
                        if der[ite] * der[ite - 1] > 0 and der[ite] * der[ite + 2] < 0 and der[ite] * der[ite + 3] < 0:
                            changeite.append(ite + 1)
                            change.append(ener[ite + 1])
    elif mode == "adapted":
        ite = 3
        while ite < len(der) - 4 and len(change) != 2:
            if len(change) == 0:
                if der[ite] > 0:
                    if 300 < ener[ite] < 1000:
                        cond1 = der[ite + 1] > 0 and der[ite + 2] > 0
                        cond2 = count[ite] < count[ite + 1] and count[ite] < count[ite + 2] and count[ite] < count[ite + 3] and count[ite] < count[ite + 4]
                        if cond1 or cond2:
                            changeite.append(ite)
                            change.append(ener[ite])
            elif len(change) == 1:
                if der[ite] < 0:
                    if 300 < ener[ite] < 1000:
                        # print(count[ite - 3], count[ite - 2], count[ite - 1], count[ite], count[ite + 1], count[ite + 2], count[ite + 3])
                        cond = count[ite] > count[ite - 1] and count[ite] > count[ite - 2] and count[ite] > count[ite - 3] and count[ite] > count[ite - 4] and count[ite] > count[ite + 1] and count[ite] > count[ite + 2] and count[ite] > count[ite + 3] and count[ite] > count[ite + 4]
                        if cond:
                            changeite.append(ite)
                            change.append(ener[ite])
            ite += 1
    else:
        raise ValueError("Use a correct mode : 'global' or 'adapted' possible")
    return change, changeite


def fit_validation(lr2, lerror, r2lim, errorlim):
    """

    :param lr2:
    :param lerror:
    :param r2lim:
    :param errorlim:
    :return:
    """
    for ite in range(len(lr2)):
        if lr2[ite] > r2lim and lerror[ite] < errorlim:
            return True
    return False


def poly_deg2(x, a, b):
    return a * x**2 + b * x

def fit_quad_calib(init_ts, final_ts, save_ratio=None, save_date=None, save_r2=None, save_error=None, save_nev=None):
    # import matplotlib.pyplot as plt
    # import numpy as np
    channels = np.array([237, 465, 599, 738, 864, 1476, 1852, 2218, 2651, 2915, 3235, 3992, 4824, 5298])
    energies = np.array([76, 186, 242, 295, 352, 609, 768, 934, 1120, 1238, 1395, 1765, 2204, 2448])
    channel_sig = np.array([44.74, 19.74, 18.73, 20.85, 23.59, 31.83, 41.65, 42.87, 50.4, 50.15, 53.22, 76.28, 192.24, 69.88])

    popt1 = curve_fit(poly_deg2, channels, energies)[0]
    # popt2 = curve_fit(poly_deg2, energies, channels)[0]
    popt2 = curve_fit(poly_deg2, energies, channels, sigma=channel_sig)[0]
    print(popt1)
    print(popt2)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    # ax.plot(channels, energies)
    ax1.errorbar(channels, energies, xerr=channel_sig)
    ax1.plot(np.linspace(0, 11000, 100), poly_deg2(np.linspace(0, 11000, 100), *popt1), color='red', label=f"{popt1[0]}x²+{popt1[1]}x")
    ax1.set(xlabel="Channel (ADC)", ylabel="Energy (keV)")
    ax1.legend()

    ax2.errorbar(energies, channels, yerr=channel_sig)
    ax2.plot(np.linspace(0, 5000, 100), poly_deg2(np.linspace(0, 5000, 100), *popt2), color='red', label=f"{popt2[0]}x²+{popt2[1]}x")
    ax2.set(xlabel="Energy (keV)", ylabel="Channel (ADC)")
    ax2.legend()
    plt.show()

    if save_ratio is not None and save_date is not None and save_r2 is not None and save_error is not None and save_nev is not None:
        smoothed3 = []
        smoothed3.append((save_ratio[0] + save_ratio[1]) / 2)
        for ite in range(1, len(save_ratio) - 1):
            smoothed3.append((save_ratio[ite - 1] + save_ratio[ite] + save_ratio[ite + 1]) / 3)
        smoothed3.append((save_ratio[-2] + save_ratio[-1]) / 2)

        no_outliers = []
        no_outliers.append(save_ratio[0])
        for ite in range(1, len(save_ratio) - 1):
            if abs(save_ratio[ite] - save_ratio[ite - 1]) > 0.09 and abs(save_ratio[ite] - save_ratio[ite + 1]) > 0.09:
                no_outliers.append((save_ratio[ite - 1] + save_ratio[ite + 1]) / 2)
            else:
                no_outliers.append(save_ratio[ite])
        no_outliers.append(save_ratio[-1])

        no_outliers_smoothed3 = []
        no_outliers_smoothed3.append((no_outliers[0] + no_outliers[1]) / 2)
        for ite in range(1, len(no_outliers) - 1):
            no_outliers_smoothed3.append((no_outliers[ite - 1] + no_outliers[ite] + no_outliers[ite + 1]) / 3)
        no_outliers_smoothed3.append((no_outliers[-2] + no_outliers[-1]) / 2)

        no_outliers_smoothed5 = []
        no_outliers_smoothed5.append((no_outliers[0] + no_outliers[1] + no_outliers[2]) / 3)
        no_outliers_smoothed5.append((no_outliers[0] + no_outliers[1] + no_outliers[2] + no_outliers[3]) / 4)
        for ite in range(2, len(no_outliers) - 2):
            no_outliers_smoothed5.append(
                (no_outliers[ite - 2] + no_outliers[ite - 1] + no_outliers[ite] + no_outliers[ite + 1] + no_outliers[ite + 2]) / 5)
        no_outliers_smoothed5.append((no_outliers[-4] + no_outliers[-3] + no_outliers[-2] + no_outliers[-1]) / 4)
        no_outliers_smoothed5.append((no_outliers[-3] + no_outliers[-2] + no_outliers[-1]) / 3)

        ite_ratio = 0
        fit_time = (np.array([int(date.split(" ")[-1]) for date in save_date]) + np.array([int(date.split(" ")[0]) for date in save_date]))/2
        interpol_ratio = []
        for ts_val in range(init_ts, final_ts+1):
            print(ts_val)
            if ts_val <= fit_time[0]:
                interpol_ratio.append(no_outliers_smoothed5[0] + (no_outliers_smoothed5[1] - no_outliers_smoothed5[0]) / (fit_time[1] - fit_time[0]) * (ts_val - fit_time[0]))
            elif ts_val >= fit_time[-1]:
                interpol_ratio.append(no_outliers_smoothed5[-2] + (no_outliers_smoothed5[-1] - no_outliers_smoothed5[-2]) / (fit_time[-1] - fit_time[-2]) * (ts_val - fit_time[-2]))
            else:
                while ite_ratio < len(fit_time) - 1 and ts_val >= fit_time[ite_ratio+1]:
                    ite_ratio += 1
                if ts_val >= fit_time[ite_ratio]:
                    interpol_ratio.append(no_outliers_smoothed5[ite_ratio] + (no_outliers_smoothed5[ite_ratio+1] - no_outliers_smoothed5[ite_ratio]) / (fit_time[ite_ratio+1] - fit_time[ite_ratio]) * (ts_val - fit_time[ite_ratio]))
                else:
                    raise IndexError("Problem during the interpolation, ts range impossible to find")


        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
        ax1.plot(fit_time, save_ratio, color="blue", label="Original")
        # ax1.scatter(fit_time, save_ratio, color="pink", s=10)
        # ax1.plot(fit_time, no_outliers, color="green")
        # ax1.plot(fit_time, smoothed3, color="orange")
        # ax1.scatter(fit_time, smoothed3, color="pink", s=10)
        ax1.plot(fit_time, no_outliers_smoothed3, color="red", label="Smoothed 3 points")
        ax1.plot(fit_time, no_outliers_smoothed5, color="orange", label="Smoothed 5 points")
        ax1.scatter(range(init_ts, final_ts+1), interpol_ratio, color="black", s=5)
        ax1.set(ylabel="511 ratio")
        ax1.legend()
        ax2.plot(fit_time, save_nev)
        ax2.set(ylabel="n events")
        ax3.plot(fit_time, save_r2)
        ax3.set(ylabel="R²")
        ax4.plot(fit_time, save_error)
        ax4.set(xlabel="Time (s)", ylabel="Relative error (%)")
        plt.show()

        with open("Kiruna_data/calib_maud/fine_calib.txt", "w") as f:
            f.write("Timestamp - Correction to be applied during the calibration to account for the temperature effect in the detector\n")
            for ite in range(len(interpol_ratio)-1):
                f.write(f"{range(init_ts, final_ts)[ite]} {interpol_ratio[ite]}\n")
            f.write(f"{final_ts} {interpol_ratio[-1]}")

# fit_quad_calib()

analysis = KirunaAnalysis()
ergcut = None
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

import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas as pd

import matplotlib as mpl

mpl.use("Qt5Agg")

def show_calib_params(calibfile):
    with open(calibfile) as f:
        data = [line.split("\t")[1:] for line in f.read().split("\n")[:-1]]

    alist = [[], [], [], [], [], [], [], []]
    chanlist = [[], [], [], [], [], [], [], []]
    names = ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7"]
    for piste, line in enumerate(data):
        for ite in range(8):
            chanlist[ite].append([float(line[ite]), piste])
            alist[ite].append(float(line[ite]))

    chanlist = np.array(chanlist)
    fig, axes = plt.subplots(2, 4)
    axs = np.array(axes).flatten()
    for itecoef, ax in enumerate(axs):
        ax.hist(alist[itecoef], bins=20, label=f"Coef {names[itecoef]}")
        ax.legend()

    fig2, axes2 = plt.subplots(2, 4)
    axs2 = np.array(axes2).flatten()
    for itecoef, ax2 in enumerate(axs2):
        ax2.scatter(chanlist[itecoef, :, 0], chanlist[itecoef, :, 1], s=5, label=f"Coef {names[itecoef]}")
        ax2.legend()
        ax2.grid(True)

calib = "./Kiruna_data/calib_dssd/nside_v2.calib"
show_calib_params(calib)


import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas as pd
from scipy.optimize import curve_fit

import matplotlib as mpl

mpl.use("Qt5Agg")

def show_dssd_temp_dependance(file_bulk="./Kiruna_data/calib_dssd/fit_bulk.txt", file_peak="./Kiruna_data/calib_dssd/fit_pedest.txt"):
    with open(file_peak) as f:
        data = f.read().split("\n")[:-1]
        peak_mean = []
        peak_mean_err = []
        for dat in data:
            peak_mean.append(dat.split("|")[0].split(" "))
            peak_mean_err.append(dat.split("|")[1].split(" "))
        peak_mean = np.array(peak_mean, dtype=np.float64)
        peak_mean_err = np.array(peak_mean_err, dtype=np.float64)
        # print(peak_mean)
        # print(peak_mean_err)

    with open(file_bulk) as f:
        data = f.read().split("\n")[:-1]
        bulk_mean = []
        bulk_mean_err = []
        for dat in data:
            bulk_mean.append(dat.split("|")[0].split(" "))
            bulk_mean_err.append(dat.split("|")[1].split(" "))
        bulk_mean = np.array(bulk_mean, dtype=np.float64)
        bulk_mean_err = np.array(bulk_mean_err, dtype=np.float64)
        # print(bulk_mean)
        # print(bulk_mean_err)

    with uproot.open("./corr_abs_dssd_file.root") as file:
        tree = file["Events"]
        df = pd.concat([chunk for chunk in tree.iterate(["time_corr_abs", "energy", "temperature"], step_size=1000000, library="pd")],
                       ignore_index=True)

    fig_temp, axes = plt.subplots(2, 3, figsize=(18, 12))
    # temps = [">10", "10-5", "5-0", "0--5", "-5--10", "<-10"]
    print(df)
    print(df[df.temperature > 10].temperature.values)
    axes[0][0].hist(df[df.temperature > 10].temperature.values, label=">10")
    axes[0][0].axvline(np.mean(df[df.temperature > 10].temperature.values), color="green")
    axes[0][1].hist(df[np.logical_and(df.temperature <= 10, df.temperature > 5)].temperature.values, label="10-5")
    axes[0][1].axvline(np.mean(df[np.logical_and(df.temperature <= 10, df.temperature > 5)].temperature.values), color="green")
    axes[0][2].hist(df[np.logical_and(df.temperature <= 5, df.temperature > 0)].temperature.values, label="5-0")
    axes[0][2].axvline(np.mean(df[np.logical_and(df.temperature <= 5, df.temperature > 0)].temperature.values), color="green")
    axes[1][0].hist(df[np.logical_and(df.temperature <= 0, df.temperature > -5)].temperature.values, label="0--5")
    axes[1][0].axvline(np.mean(df[np.logical_and(df.temperature <= 0, df.temperature > -5)].temperature.values), color="green")
    axes[1][1].hist(df[np.logical_and(df.temperature <= -5, df.temperature > -10)].temperature.values, label="-5--10")
    axes[1][1].axvline(np.mean(df[np.logical_and(df.temperature <= -5, df.temperature > -10)].temperature.values), color="green")
    axes[1][2].hist(df[df.temperature <= -10].temperature.values, label="<-10")
    axes[1][2].axvline(np.mean(df[df.temperature <= -10].temperature.values), color="green")
    for axs in axes:
        for ax in axs:
            ax.set(xlabel="Temperature", ylabel="Number of event")
            ax.legend()
    plt.show()

    xs = [np.mean(df[df.temperature > 10].temperature.values),
          np.mean(df[np.logical_and(df.temperature <= 10, df.temperature > 5)].temperature.values),
          np.mean(df[np.logical_and(df.temperature <= 5, df.temperature > 0)].temperature.values),
          np.mean(df[np.logical_and(df.temperature <= 0, df.temperature > -5)].temperature.values),
          np.mean(df[np.logical_and(df.temperature <= -5, df.temperature > -10)].temperature.values),
          np.mean(df[df.temperature <= -10].temperature.values)]

    def affine(x, a, b):
        return a*x+b

    apeak = []
    bpeak = []
    abulk = []
    bbulk = []

    fig_fitpeak, axes1 = plt.subplots(4, 8, figsize=(18, 18))
    for ite in range(32):
        popt, pcov = curve_fit(affine, xs, peak_mean[ite], sigma=peak_mean_err[ite])[:2]
        apeak.append(popt[0])
        bpeak.append(popt[1])
        axes1[int(ite/8)][ite%8].errorbar(xs, peak_mean[ite], yerr=peak_mean_err[ite], color="blue")
        axes1[int(ite/8)][ite%8].plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="red")
        axes1[int(ite/8)][ite%8].set(title=f"{round(popt[0], 2)} x + {round(popt[1], 2)}", xlabel="Temperature", ylabel="Peak position(ADC)")
    fig_fitpeak.set_constrained_layout(True)

    plt.show()

    fig_fitbulk, axes2 = plt.subplots(4, 8, figsize=(18, 18))
    for ite in range(32):
        popt, pcov = curve_fit(affine, xs, bulk_mean[ite], sigma=bulk_mean_err[ite])[:2]
        abulk.append(popt[0])
        bbulk.append(popt[1])
        axes2[int(ite/8)][ite%8].errorbar(xs, bulk_mean[ite], yerr=bulk_mean_err[ite], color="red")
        axes2[int(ite/8)][ite%8].plot(np.linspace(-15, 15, 100), affine(np.linspace(-15, 15, 100), *popt), color="blue")
        axes2[int(ite/8)][ite%8].set(title=f"{round(popt[0], 2)} x + {round(popt[1], 2)}", xlabel="Temperature", ylabel="Bulk position(ADC)")
    fig_fitbulk.set_constrained_layout(True)

    plt.show()

    # with open("./Kiruna_data/calib_dssd/dssd_temp_corr_coefs", "w") as f:
    #     f.write("std::vector<Double_t> peak_a = {")
    #     for value in apeak[:-1]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(apeak[-1], 6)}")
    #     f.write("};\n")
    #     f.write("std::vector<Double_t> peak_b = {")
    #     for value in bpeak[:-1]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bpeak[-1], 6)}")
    #     f.write("};\n")
    #     f.write("std::vector<Double_t> bulk_a = {")
    #     for value in abulk[:-1]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(abulk[-1], 6)}")
    #     f.write("};\n")
    #     f.write("std::vector<Double_t> bulk_b = {")
    #     for value in bbulk[:-1]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bbulk[-1], 6)}")
    #     f.write("};\n\n")
    #
    #     f.write("std::vector<std::vector<Double_t>> peak_adc = {\n")
    #     f.write("{")
    #     for value in peak_mean[:-1, 0]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(peak_mean[-1, 0], 6)}")
    #     f.write("},\n{")
    #     for value in peak_mean[:-1, 1]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(peak_mean[-1, 1], 6)}")
    #     f.write("},\n{")
    #     for value in peak_mean[:-1, 2]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(peak_mean[-1, 2], 6)}")
    #     f.write("},\n{")
    #     for value in peak_mean[:-1, 3]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(peak_mean[-1, 3], 6)}")
    #     f.write("},\n{")
    #     for value in peak_mean[:-1, 4]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(peak_mean[-1, 4], 6)}")
    #     f.write("},\n{")
    #     for value in peak_mean[:-1, 5]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(peak_mean[-1, 5], 6)}")
    #     f.write("}\n};\n\n")
    #
    #     f.write("std::vector<std::vector<Double_t>> bulk_adc = {\n")
    #     f.write("{")
    #     for value in bulk_mean[:-1, 0]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bulk_mean[-1, 0], 6)}")
    #     f.write("},\n{")
    #     for value in bulk_mean[:-1, 1]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bulk_mean[-1, 1], 6)}")
    #     f.write("},\n{")
    #     for value in bulk_mean[:-1, 2]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bulk_mean[-1, 2], 6)}")
    #     f.write("},\n{")
    #     for value in bulk_mean[:-1, 3]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bulk_mean[-1, 3], 6)}")
    #     f.write("},\n{")
    #     for value in bulk_mean[:-1, 4]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bulk_mean[-1, 4], 6)}")
    #     f.write("},\n{")
    #     for value in bulk_mean[:-1, 5]:
    #         f.write(f"{round(value, 6)}, ")
    #     f.write(f"{round(bulk_mean[-1, 5], 6)}")
    #     f.write("}\n};")

show_dssd_temp_dependance()
show_dssd_temp_dependance(file_bulk="./Kiruna_data/calib_dssd/fit_corr_bulk.txt", file_peak="./Kiruna_data/calib_dssd/fit_corr_pedest.txt")