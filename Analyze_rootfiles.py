import uproot
import pandas as pd
import numpy as np
from scipy.fftpack import rfft, rfftfreq
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt
import matplotlib as mpl
from memory_profiler import memory_usage
from datetime import datetime
mpl.use('TkAgg')


print(f"Memory at first : {memory_usage(-1)} MB")

# Mean flux of 0.3 ph/s but emission every 0.03 sec : mean num of ph emited per pulsation : 0.3*0.03 = 0.009
crab_flux = 0.003
crab_period = 33.9e-3
crab_frequency = round(1/33.9e-3, 3)

def format_axis(func_fftvals, len_signal, bin_timescale):
    x_axis = rfftfreq(len_signal, bin_timescale)
    y_axis = np.abs(func_fftvals)
    index = np.where(x_axis > 0.1, True, False)
    return x_axis[index], y_axis[index]


def crab_draw(xlist, mean, bin_timescale):
    added_signal = []
    for centroid in xlist:
        if centroid < bin_timescale/2:
            added_signal.append(poisson.rvs(mean))
        elif centroid > crab_period-bin_timescale/2:
            added_signal.append(poisson.rvs(mean))
        else:
            added_signal.append(0)
    added_signal = np.array(added_signal)
    return added_signal


def run_fft(data_signal, fft_title, namefig=None, bin_timescale=1e-2, vline=None):
    print(f" Running a Fourier Fast Transform analysis : {fft_title}")
    # Calcul de la TF rapide (FFT)
    func_fftvals = rfft(data_signal)
    # print(f"Memory in the fft after applying the fft : {memory_usage(-1)} MB")

    new_x, new_y = format_axis(func_fftvals, len(data_signal), bin_timescale)
    # print(f"Memory in the fft after making the axis for the plots : {memory_usage(-1)} MB")

    title = f"{fft_title}, bining timescale : {bin_timescale}"
    # Visualisation de la TF
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    if vline is not None:
        ax.axvline(vline, color="red", label="Crab pulsar frequency")
    ax.plot(new_x, new_y, label="FFT of the signal")
    ax.set(xlabel='Frequency [Hz]', ylabel='Amplitude', title=title)
    ax.legend()
    if namefig is not None:
        plt.savefig(namefig)
    plt.show()


def get_bined_times(datafile, fft_title, crab_fl, bin_timescale=1e-2, add_signal=False, emptyness=None):
    with uproot.open(datafile) as file:
        # file is there a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        # Maud tree keys :
        # print(maud_tree.keys())

        # convert the tree in dataframe
        df = tree.arrays(library="pd")
        print(f"Memory after opening rootfile as a df ({len(df)} events): {memory_usage(-1)} MB")
        # print("n events : ", len(df))

        # print(f"Memory after converting the root file in df : {memory_usage(-1)} MB")

        # Afficher le DataFrame
        # print(df_maud.head())
        time_array = df.time_corr_abs.values

        # tinit, tfinal = time_array[0] + 1719079200, time_array[-1] + 1719079200
        # print(tinit, tfinal)
        # print(int(tinit), int(tfinal))
        # print(datetime.utcfromtimestamp(int(tinit)), datetime.utcfromtimestamp(int(tfinal)))
        #
        # stop
        print("===============================================================================================================")
        print(f" Preparing the fft : {fft_title}")
        init_time = round(int(time_array[0]/bin_timescale) * bin_timescale, 3)
        end_time = round((int(time_array[-1]/bin_timescale) + 1) * bin_timescale, 3)
        n_val = int((end_time - init_time) / bin_timescale + 1)
        bins = np.linspace(init_time, end_time, n_val)
        print(f"  == init time      : {init_time}")
        print(f"  == end_time       : {end_time}")
        print(f"  == number of bins : {n_val - 1}")
        # print(bins)
        # print(len(bins))
        # print(n_val)
        data_signal = np.histogram(time_array, bins=bins)[0]
        # print(len(data_signal))
        if add_signal:
            centroids = np.around(bins[:-1] % crab_period, 5)
            added_signal = crab_draw(centroids, crab_fl, bin_timescale)
            if emptyness is not None:
                correction_list = []
                for ite in range(len(bins)-1):
                    verif_bool = False
                    for inf_lim, sup_lim in emptyness:
                        if bins[ite] >= inf_lim and bins[ite+1] <= sup_lim:
                            verif_bool = True
                    if verif_bool:
                        correction_list.append(0)
                    else:
                        correction_list.append(1)
                print(f"before correction : {np.sum(added_signal)}")
                added_signal = added_signal * np.array(correction_list)
                print(f"after correction : {np.sum(added_signal)}")

            # print(f"nb counts add sig = {np.sum(added_signal[30079714:30379714])}")
            # print(f"Expected count = 900")
            # print(added_signal[:40])
            # print(sig2[:40])
            # print(data_signal[:40])
            data_signal = data_signal + added_signal
        # print(data_signal[:40])

        # print(n_val, len(data_signal))
        print(f"Memory in the fft after creating the histograms (before releasing memory) : {memory_usage(-1)} MB")
        return data_signal


def get_bined_selected_times(datafile, fft_title, time_ranges, crab_fl, bin_timescale=1e-2, add_signal=False, emptyness=None):
    if type(time_ranges) is not list and type(time_ranges[0]) is not list:
        raise TypeError(f"time_ranges must be a list of list : {type(time_ranges), type(time_ranges[0])}")
    with uproot.open(datafile) as file:
        tree = file["Events"]
        # Maud tree keys :
        # print(f"Memory after opening rootfile : {memory_usage(-1)} MB")

        # convert the tree in dataframe
        df = tree.arrays(library="pd")
        print(f"Memory after opening rootfile as a df ({len(df)} events): {memory_usage(-1)} MB")
        # print("n events : ", len(df))

        time_array = df.time_corr_abs.values

        # for vals in time_ranges:
        #     print(datetime.utcfromtimestamp(vals[0]), datetime.utcfromtimestamp(vals[1]))
        #
        # stop
        print("===============================================================================================================")
        print(f" Preparing the fft : {fft_title}")
        init_time = round(int(time_array[0]/bin_timescale) * bin_timescale, 3)
        end_time = round((int(time_array[-1]/bin_timescale) + 1) * bin_timescale, 3)
        n_val = int((end_time - init_time) / bin_timescale + 1)
        bins = np.linspace(init_time, end_time, n_val)
        print(f"  == init time      : {init_time}")
        print(f"  == end_time       : {end_time}")
        print(f"  == number of bins : {n_val - 1}")

        visible_array = []
        for val in time_array:
            seen = False
            for time_range in time_ranges:
                if time_range[0] <= val + 1719079200 <= time_range[1]:
                    seen = True
                # else:
                #     print(val + 1719079200 - time_range[0], time_range[1] - (val + 1719079200))
            if seen:
                visible_array.append(val)
        visible_array = np.array(visible_array)
        # print(visible_array)
        data_signal = np.histogram(visible_array, bins=bins)[0]
        if add_signal:
            centroids = np.around(bins[:-1] % crab_period, 5)
            # Mean flux of 0.3 ph/s but emission every 0.03 sec : mean num of ph emited per pulsation : 0.3*0.03 = 0.009
            added_signal = crab_draw(centroids, crab_fl, bin_timescale)
            if emptyness is not None:
                correction_list = []
                for ite in range(len(bins)-1):
                    verif_bool = False
                    for inf_lim, sup_lim in emptyness:
                        if bins[ite] >= inf_lim and bins[ite+1] <= sup_lim:
                            verif_bool = True
                    if verif_bool:
                        correction_list.append(0)
                    else:
                        correction_list.append(1)
                print(f"before correction : {np.sum(added_signal)}")
                added_signal = added_signal * np.array(correction_list)
                print(f"after correction : {np.sum(added_signal)}")

            data_signal = data_signal + added_signal
        return data_signal


def get_fake_signal(fft_title, crab_fl, bin_timescale=1e-2, add_signal=False, emptyness=None):
    print("===============================================================================================================")
    print(f" Preparing the fft : {fft_title}")
    # == Full data set
    init_time = 19202.86
    end_time = 323335.52
    # == Short data set
    # init_time = 95332.18
    # end_time = 104507.18
    n_val = int((end_time - init_time) / bin_timescale + 1)
    bins = np.linspace(init_time, end_time, n_val)
    print(f"  == init time      : {init_time}")
    print(f"  == end_time       : {end_time}")
    print(f"  == number of bins : {n_val - 1}")

    centroids = np.around(bins[:-1] % crab_period, 5)
    # Total flux between 320000-323000 s : 506449
    # Source flux of 0.3 ph/s so 900 ph over this duration
    # >>> background flux integrated over 3000s : 505549  >>>> mean flux of 168.52 ph/s but the bin is not 1s wide so timescale corrects
    bkgflux = 168.52 * bin_timescale
    bkg_signal = np.array([poisson.rvs(bkgflux) for ite in range(len(bins)-1)])
    if add_signal:
        # Mean flux of 0.3 ph/s but emission every 0.03 sec : mean num of ph emited per pulsation : 0.3*0.03 = 0.009
        added_signal = crab_draw(centroids, crab_fl, bin_timescale)
        bkg_signal = bkg_signal + added_signal
    if emptyness is not None:
        correction_list = []
        for ite in range(len(bins) - 1):
            verif_bool = False
            for inf_lim, sup_lim in emptyness:
                if bins[ite] >= inf_lim and bins[ite + 1] <= sup_lim:
                    verif_bool = True
            if verif_bool:
                correction_list.append(0)
            else:
                correction_list.append(1)
        print(f"before correction : {np.sum(bkg_signal)}")
        bkg_signal = bkg_signal * np.array(correction_list)
        print(f"after correction : {np.sum(bkg_signal)}")

    # print(n_val, len(data_signal))
    print(f"Memory in the fft after creating the histograms (before releasing memory) : {memory_usage(-1)} MB")
    return bkg_signal

# print(f"Memory after creating the functions : {memory_usage(-1)} MB")

# EXAMPLE
makeex = False
if makeex:
    timescale = 1e-2
    duration = 10
    numval = int(duration / timescale)
    tmps = np.linspace(0, duration, numval, endpoint=False)
    sine = np.sin(tmps*2*np.pi*2) + np.sin(tmps*2*np.pi*1) + np.sin(tmps*2*np.pi*29.5)

    # plt.plot(tmps, sine)
    # plt.show()

    # Calcul de la TF rapide (FFT)
    fftvals = rfft(sine)

    # Visualisation de la TF
    plt.plot(rfftfreq(numval, timescale), np.abs(fftvals))
    plt.xlabel('Fréquence')
    plt.ylabel('Amplitude')
    plt.show()

    plt.show()

makeex2 = False
if makeex2:
    def entry(time):
        gap = min(time%crab_period, crab_period - time%crab_period)
        return poisson.rvs(norm.pdf(gap, loc=0, scale=0.008))

    timescale = 1e-2
    duration = 1
    numval = int(duration / timescale)
    tmps = np.linspace(0, duration, numval, endpoint=False)
    signal = np.array([entry(tmp) for tmp in tmps])
    ref_sig = 25 + 25 * np.cos(tmps*2*np.pi/crab_period)


    # plt.plot(tmps, signal)
    # plt.plot(tmps, ref_sig, color='red')
    # plt.show()

    # Calcul de la TF rapide (FFT)
    fftvals = rfft(signal)

    # Visualisation de la TF
    plt.plot(rfftfreq(numval, timescale), np.abs(fftvals))
    plt.xlabel('Fréquence')
    plt.ylabel('Amplitude')
    plt.show()

############################################################################################################
# For information on uproot :
# https://uproot.readthedocs.io/en/latest/basic.html
############################################################################################################
check_gaps_in_data = False
if check_gaps_in_data:
    with uproot.open("corr_abs_maud_file.root") as file:
        # file is there a TDirectory, you can display the content with :
        # print(file.keys())
        # You can have even more detail on what the key is with :
        # print(file.classnames())
        # Accessing the tree
        tree = file["Events"]
        # Maud tree keys :
        # print(maud_tree.keys())

    # convert the tree in dataframe
    df = tree.arrays(library="pd")
    print("n events : ", len(df))

    time_array = df.time_corr_abs.values
    ite_list = []
    for ite in range(len(time_array) - 1):
        if time_array[ite + 1] - time_array[ite] > 0.5:
            ite_list.append([time_array[ite], time_array[ite + 1]])
    print(ite_list)
    check_times = False
    if check_times:
        # Afficher le DataFrame
        # print(df_maud.head())
        gps_from_pps_array = df.pps_cpt_corr_abs.values + 1719079200
        diff_list = []
        for ite in range(len(time_array) - 1):
            if time_array[ite + 1] - time_array[ite] > 1:
                print(f"Diff : {time_array[ite + 1] - time_array[ite]}   ite {ite}-{ite + 1}/{len(time_array)}")
                diff_list.append(ite)
                date_time_local = datetime.fromtimestamp(gps_from_pps_array[ite])
                date_time_localp1 = datetime.fromtimestamp(gps_from_pps_array[ite+1])
                print(f"Date et heure locale : {date_time_local} - {date_time_localp1}\n")
            if time_array[ite + 1] <= time_array[ite]:
                print("ERROR")

        print(f"Num of vals in a row : {diff_list[0]}")
        for ite in range(len(diff_list)-1):
            print(f"Num of vals in a row : {diff_list[ite+1]-diff_list[ite]}")
        print(f"Num of vals in a row : {len(time_array)-diff_list[-1]}")

times_of_emptyness = [[19237.272103425, 21130.70923795], [22409.7524636, 22434.4893218], [26444.113075475, 28864.92999595],
                      [64440.98104165, 64441.6488771], [85554.701873725, 85891.042379775], [248483.597786375, 248630.742785725],
                      [248800.05419025, 248883.775494775], [249157.926650425, 249382.772406325], [251492.077928625, 251649.213788825],
                      [260324.563231225, 260327.554940825], [261335.77597295, 261336.297486475], [323258.856362225, 323304.646581125]]

# which_sim = None
# which_sim = "root"
# which_sim = "improvedroot"
which_sim = "random"
maudfile = "corr_abs_maud_file.root"
# maudfile = "corr_abs_maud_file_short.root"
timescale = 1e-2
if which_sim == "root":
    # print(f"Memory before opening rootfile : {memory_usage(-1)} MB")
    # Opening ROOT file
    print(f"Memory before getting the signal : {memory_usage(-1)} MB")
    name = "FFT with corrected time"
    signal = get_bined_times(maudfile, name, crab_flux, bin_timescale=timescale, add_signal=False, emptyness=times_of_emptyness)
    # signal = get_bined_times(maudfile, "FFT with corrected time", bin_timescale=timescale, add_signal=True, emptyness=None)
    # signal = get_bined_times(maudfile, "FFT with corrected time", bin_timescale=timescale, add_signal=True, emptyness=times_of_emptyness)
    print(f"Memory after getting the signal, before the fft : {memory_usage(-1)} MB")
    run_fft(signal, name, bin_timescale=timescale, vline=crab_frequency)
    # run_fft(uncorr_time_array, "FFT with uncorrected time")
    # run_fft(uncorr_time_array_abs, "FFT with uncorrected time with common time origin")
    print(f"Memory after the fft : {memory_usage(-1)} MB")
elif which_sim == "random":
    print(f"Memory before getting the signal : {memory_usage(-1)} MB")
    name = "FFT with fake signal"
    # signal = get_fake_signal("FFT with fake signal", bin_timescale=timescale, add_signal=False, emptyness=None)
    # signal = get_fake_signal("FFT with fake signal", bin_timescale=timescale, add_signal=False, emptyness=times_of_emptyness)
    # signal = get_fake_signal("FFT with fake signal", bin_timescale=timescale, add_signal=True, emptyness=None)
    signal = get_fake_signal(name, crab_flux, bin_timescale=timescale, add_signal=True, emptyness=times_of_emptyness)
    print(f"Memory after getting the signal, before the fft : {memory_usage(-1)} MB")
    run_fft(signal, name, bin_timescale=timescale, vline=crab_frequency)
    print(f"Memory after the fft : {memory_usage(-1)} MB")
if which_sim == "improvedroot":
    # Opening ROOT file
    print(f"Memory before getting the signal : {memory_usage(-1)} MB")
    # alt > 5°
    # kept_times = [["23/2h", "23/20h30"], ["24/2h40", "24/23h40"], ["25/3h50", "26/2h"], ["2h45", "fin"]]
    # kept_times, fignum, deg_lim = [[1719108000, 1719174600], [1719196800, 1719272400], [1719287400, 1719367200], [1719377100, 1719446400]], "all", 5
    # kept_times, fignum, deg_lim = [[1719108000, 1719174600]], 1, 5
    # kept_times, fignum, deg_lim = [[1719196800, 1719272400]], 2, 5
    # kept_times, fignum, deg_lim = [[1719287400, 1719367200]], 3, 5
    # kept_times, fignum, deg_lim = [[1719377100, 1719446400]], 4, 5
    # # alt > 15°
    # # kept_times = [["23/4h20", "23/18h10"], ["24/5h", "24/21h"], ["25/6h40", "25/22h50"], ["26/8h", "fin"]]
    # kept_times, fignum, deg_lim = [[1719116400, 1719166200], [1719205200, 1719262800], [1719297600, 1719355800], [1719388800, 1719446400]], "all", 15
    # kept_times, fignum, deg_lim = [[1719116400, 1719166200]], 1, 15
    # kept_times, fignum, deg_lim = [[1719205200, 1719262800]], 2, 15
    # kept_times, fignum, deg_lim = [[1719297600, 1719355800]], 3, 15
    kept_times, fignum, deg_lim = [[1719388800, 1719446400]], 4, 15
    name = "FFT with corrected time"
    signal = get_bined_selected_times(maudfile, name, kept_times, crab_flux, bin_timescale=timescale, add_signal=False, emptyness=None)
    print(f"Memory after getting the signal, before the fft : {memory_usage(-1)} MB")
    run_fft(signal, name, namefig=f"window_{fignum}_{deg_lim}°", bin_timescale=timescale, vline=crab_frequency)
    print(f"Memory after the fft : {memory_usage(-1)} MB")

