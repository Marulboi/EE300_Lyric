import numpy as np

import matplotlib.pyplot as plt

import scipy.signal as signal

from scipy.signal import savgol_filter, find_peaks
import pandas as pd

from scipy.interpolate import interp1d

from scipy.io import loadmat ,savemat
import glob , os

sampling_rate = 2048.0     # Sampling rate in Hz
prominenceForPeaks = 0.1
highttreshold = 0.7
msTime= 50

def find_csv_mat_pairs(base_path):
    signals_csv_folder = os.path.join(base_path, "SignalsCSV")
    signals_beats_folder = os.path.join(base_path, "Signals_Beats")

    # Find all CSV files
    csv_files = glob.glob(os.path.join(signals_csv_folder, "*.csv"))

    pairs = []

    for csv_path in csv_files:
        # Get base name without extension
        csv_name = os.path.splitext(os.path.basename(csv_path))[0]
        
        # Expected .mat file name
        mat_name = f"{csv_name}_Signal_select_Beat.mat"
        mat_path = os.path.join(signals_beats_folder, mat_name)

        if os.path.exists(mat_path):
            pairs.append((csv_path, mat_path))
        else:
            print(f"[!] Missing .mat file for: {csv_name}")

    return pairs





# Load the filtered signal from CSV and transpose it

base_folder = r"C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig5"
pairs = find_csv_mat_pairs(base_folder)


for filename, matfilename in pairs:
    print(f"CSV: {os.path.basename(filename)}")
    print(f"MAT: {os.path.basename(matfilename)}\n")
    try:
    
        mat_data = loadmat(matfilename)
        mat_data = mat_data['Signal_select_Beat']

        QRSonset = mat_data[0][0][0] #Majik(magic) numbers I know but it is to acces the data which is alyas structered this way
        QRSonset = QRSonset.flatten()
        QRSrange = mat_data[0][0][2][0][0]


        df = pd.read_csv(filename, sep=',')
        df = df.T  # Rows become channels, columns become time samples

        # List of bad leads (channels to exclude from analysis)
        #Badlead_Ve = [1,15,	16,	25,	46,	53,	71,	72,	74,	78,	84,	96,	97, 105, 152, 153,]

        # Remove bad leads from the dataframe
        #df.drop(Badlead_Ve, axis=1, inplace=True)

        # Create a time vector assuming a sampling rate of 2048 Hz
        time = np.linspace(0, df[0].size / 2048, df[0].size)

        min_distance_qrs_ms = 300  # Minimum time between two QRS complexes (in ms)
        sampling_rate = 2048.0     # Sampling rate in Hz

        peak_counts = []
        ecgleads = []




        # Loop over all leads (columns)
        for lead in df.columns:
            y = df[lead][:].copy().to_numpy()  # Get the signal
            
            sampling_rate = sampling_rate
            # Detect peaks
            peaks, _ = find_peaks(
            y, 
            prominence=prominenceForPeaks, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=highttreshold * np.max(y)  # absouute minthreshold
            )

            troughs, _ = find_peaks(
            -y, 
            prominence=prominenceForPeaks, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=highttreshold * np.max(-y)  # absolute minthreshold
            )

            # Store number of detected peaks
            peak_counts.append(len(peaks)+len(troughs))

            # Find the index of the lead with median number of peaks
        median_index = np.argsort(peak_counts)[round((len(peak_counts)*7)/8)]
        
        n_ecg = df.columns[median_index]



        x = df[n_ecg][(QRSonset[0]):(QRSonset[0] + QRSrange)-1].copy().to_numpy()   # Make a copy of part of the signal for processing 
        time = np.linspace(0, x.size / 2048, x.size)
        # Define parameters for peak detection
        min_distance_qrs_ms = 300  # Minimum time between two QRS complexes (in ms)
        sampling_rate = 2048.0     # Sampling rate in Hz

        # Detect negative peaks (likely R-troughs in ECG) with minimum prominence and spacing
        time_low_local_peaks, _ = find_peaks(
            -x, 
            prominence=prominenceForPeaks, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=0.3 * np.max(-x)  
        )

        # Detect positive peaks (likely R-peaks in ECG)
        time_high_local_peaks, _ = find_peaks(
            x, 
            prominence=prominenceForPeaks, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=0.3 * np.max(x) 
        )

        # plt.plot(time, x)  # Raw signal
        # plt.plot(time[time_low_local_peaks], x[time_low_local_peaks], "o")  # Troughs
        # plt.plot(time[time_high_local_peaks], x[time_high_local_peaks], "x")  # Peaks
        # plt.show()
        # print(time_high_local_peaks)

        if (len(time_high_local_peaks)>len(time_low_local_peaks)):
            peakoffset = time_low_local_peaks[0]+10
            newQRSOnset = troughs - peakoffset
        else:
            peakoffset = time_high_local_peaks[0]+10
            newQRSOnset = peaks - peakoffset

        newQRSOnset = newQRSOnset[1:-1]

        # plt.plot(time, x)  # Raw signal
        # plt.plot(time[newQRSOnset], x[newQRSOnset], "o")  # QRSonsetnew
        # plt.plot(time[QRSonset], x[QRSonset], "^")  # Troughs
        # plt.plot(time[newQRSOnset+QRSrange], x[newQRSOnset+QRSrange], "*")  # Troughs
        # plt.plot(time[peaks], x[peaks], "x")  # Peaks
        # plt.show()
        
        # === Prepare new save path ===
        base_folder = os.path.dirname(os.path.dirname(matfilename))  # Go up from Signals_Beats
        save_dir = os.path.join(base_folder, "Signal_select_auto")
        os.makedirs(save_dir, exist_ok=True)

        # === Prepare new struct with appended field ===
        # Read original struct again (already done earlier as mat_data)
        original_struct = mat_data[0][0]  # This is a numpy.void struct

        # Convert the original MATLAB struct to a Python dictionary
        struct_fields = {
            'QRS_Onset': original_struct['QRS_Onset'].astype('double'),
            'Beat_range': original_struct['Beat_range'].astype('double'),
            'QRS_range': original_struct['QRS_range'].astype('double'),
            'T_Onset': original_struct['T_Onset'].astype('double'),
            'BeatOffset': original_struct['BeatOffset'].astype('double'),
            'NoiseRange': original_struct['NoiseRange'].astype('double'),
            'QRS_Onset_auto': newQRSOnset.reshape(-1, 1).astype('double') 
        }

        # Save to the new .mat file
        new_matfile_path = os.path.join(save_dir, os.path.basename(matfilename))
        savemat(new_matfile_path, {'Signal_select_Beat': struct_fields})

    except Exception as e:
        print(f"Skipping {filename} due to error: {e}")
        continue
