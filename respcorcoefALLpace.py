
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter, find_peaks
from scipy.interpolate import interp1d
from scipy.signal import find_peaks                                              
import os

def norm_df(df):
    return (df - df.min()) / (df.max() - df.min())

def corr_with_pca(series, pca_series):
    # Ensure both are pandas Series and aligned
    pca_series = pd.Series(pca_series, index=series.index)
    return series.corr(pca_series)



folder_path = r"C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\SignalsCSV"

datafolder = os.path.split(folder_path)[0]
output_folder = os.path.join(folder_path,'Corralationswithlambda')
os.makedirs(output_folder, exist_ok=True)  # creates folder if it doesn't exist
output_folder = os.path.join(folder_path,'Corralationswithlambda')
#lead selector parameters
cutoffbottom = 0.4
cutofftop = 0.2

PigCSVDataPathList = []

for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
    if os.path.isfile(file_path):
        PigCSVDataPathList.append(file_path)
print(PigCSVDataPathList)



for filepath in PigCSVDataPathList:
    try:
        # Path to filtered respiration data files and list of bad channels
        os.path.splitext(os.path.basename(filepath))[0]



        temp = os.path.join(os.path.split(folder_path)[0],'Signal_select_auto','respcorcof',(os.path.splitext(os.path.basename(filepath))[0]+'_Respcheck.csv'))

        dfdata = pd.read_csv(temp, sep=',')

        temp= ''
        dfdata['data2'][1:]




        # Load the filtered signal from CSV and transpose it
        df = pd.read_csv(filepath, sep=',')
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
            

            # Detect peaks
            peaks, _ = find_peaks(
            y, 
            prominence=0.1, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=0.4 * np.max(y)  # absouute minthreshold
            )

            troughs, _ = find_peaks(
            -y, 
            prominence=0.1, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=0.4 * np.max(y)  # absouute minthreshold
            )

            if True :#(len(peaks)>2 and len(troughs)>2):
            # Store number of detected peaks
                peak_counts.append(len(peaks))
                ecgleads.append(y)
        # Find the index of the lead with median number of peaks
        median_index = np.argsort(peak_counts)[(len(peak_counts)*2// 3)]
        median_ecg_signal = ecgleads[median_index]

            


        n_ecg = df.columns[median_index]
        print("---------------------------------------------")
        print(f"File: {filepath} ECG lead selceted {n_ecg}")

        # Select channel number 2 (e.g., one ECG lead) from the dataframe

        x = df[n_ecg][:].copy()  # Make a copy of the signal for processing

        # Define parameters for peak detection

        # Detect negative peaks (likely R-troughs in ECG) with minimum prominence and spacing
        time_low_peaks, _ = find_peaks(
            -x, 
            prominence=0.1, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=0.4 * np.max(x)  # seuil minimum en valeur absolue
        )

        # Detect positive peaks (likely R-peaks in ECG)
        time_high_peaks, _ = find_peaks(
            x, 
            prominence=0.1, 
            distance=min_distance_qrs_ms / 1000 * sampling_rate,
            height=0.4 * np.max(x)  # seuil minimum en valeur absolue
        )

        #Ensure the signal starts with a low peak (trough)
        if time_low_peaks[0] > time_high_peaks[0]:  
            time_high_peaks = time_high_peaks[1:]  # Remove first high peak if it occurs before the first trough
            print("first high_peaks removed")

        # Ensure the signal ends with a high peak
        if time_low_peaks[-1] > time_high_peaks[-1]:
            print("last low_peaks removed")
            time_low_peaks = time_low_peaks[:-1]  # Remove last low peak if it occurs after the last high peak

        # plt.savefig(os.path.join( folder_path,'png','ECG', os.path.basename(filepath).split('.')[0]+'.png'), dpi=300, bbox_inches='tight')
        # Find and select leads for ECGI with the desired variances
        varlist = []
        for i in df.columns.tolist():
            n_ecg = i
            x = df[n_ecg][:].copy()  # Make a copy of the signal for processing

            # Define parameters for peak detection
            min_distance_qrs_ms = 300  # Minimum time between two QRS complexes (in ms)
            sampling_rate = 2048.0     # Sampling rate in Hz
            # Detect positive peaks (likely R-peaks in ECG)
            time_high_peaks2, _ = find_peaks(
                x, 
                prominence=0.1, 
                distance=min_distance_qrs_ms / 1000 * sampling_rate,
                height=0.4 * np.max(x)  # seuil minimum en valeur absolue
            )
            varlist.append((np.var(x.iloc[time_high_peaks2]),i))

        vararray = np.array(varlist)
        sorted_arr = vararray[np.argsort(vararray[:, 0])]
        n = len(sorted_arr)
        bot_index = int(n * (cutoffbottom))  
        top_index = int(n * (cutofftop))
        top_half = sorted_arr[top_index:bot_index]
        top_half = np.sort(top_half[:, 1].astype(int))
        print(top_half)

        selected_columns = df[top_half]


        # Define the number of points corresponding to a 60ms window
        # 50 samples at 2048 Hz ≈ 50 milliseconds
        nPoint50ms = int(50 * 2048 / 1000)

        # Assign the signal dataframe to a new variable
        signals = selected_columns

        # Ensure that rows are channels and columns are time samples
        if signals.shape[1] < signals.shape[0]:
            signals = signals.T

        # Convert the DataFrame to a NumPy array of floats for numerical operations
        signals = signals.to_numpy(dtype=float)

        # Initialize matrix to store RR amplitudes:
        # Shape = (60 channels × number of R-peaks)
        RRAmplitude = np.zeros((signals.shape[0], len(time_high_peaks)))

        # Initialize a placeholder vector for PCA projection (3D)
        vector = np.zeros((3, len(time_high_peaks)))

        # Restrict the signal matrix to the first 60 channels

        # Loop through each detected R-peak (time_high_peaks)
        for iPeak, localPeak in enumerate(time_high_peaks):
            
            # For each channel, get the max absolute value in a ±60ms window around the peak
            localMax = [
                np.abs(sig[max(localPeak - nPoint50ms,0) : min(localPeak + nPoint50ms,len(sig))]).max() 
                for sig in signals
            ]

            # Store the amplitudes in the corresponding column
            RRAmplitude[:, iPeak] = localMax[0:(signals.shape[0])]

        # Create a processing pipeline:
        # 1. Standardize features across channels
        # 2. Apply PCA to reduce to 3 components
        pca = Pipeline([
            ('Scaler', StandardScaler()),
            ('PCA', PCA(n_components=3))
        ])

        # Fit PCA on the transposed RRAmplitude (features in columns) and transform the data
        pca_RRamp = pca.fit_transform(RRAmplitude.T)


        # Select which PCA component to use (e.g., 0 = first principal component)
        n_component = 0

        # Normalize the selected PCA component to the range [-1, 1]
        scaler = MinMaxScaler(feature_range=(-1, 1))
        pca_RRamp_normalized = scaler.fit_transform(
            pca_RRamp[:, n_component].reshape(-1, 1)  # Reshape to 2D for scaler
        ).flatten()  # Flatten back to 1D after scaling

        """ # Create a large figure and plotting axis
        ax = plt.figure(figsize=[20, 10]).subplots()

        # Plot the normalized PCA projection against the time of R-peaks
        ax.plot(time[time_high_peaks], pca_RRamp_normalized)

        # Set Y-axis limits to match the normalization range
        ax.set_ylim([-1, 1])

        """

        reference_track = savgol_filter(pca_RRamp[:,n_component], window_length=6, polyorder=2)
        # plt.plot(reference_track)


        # Create time vector corresponding to PCA values (at each R-peak)
        time_reference = time[time_high_peaks]

        # Create an interpolation function to reconstruct a continuous signal
        interp_func = interp1d(
            time_reference,         # Known time points (where PCA was computed)
            reference_track,        # Corresponding PCA values
            kind='linear',          # Interpolation type: 'linear' 
            fill_value="extrapolate"  # Allow extrapolation beyond the last known time point
        )

        # Use the original time vector at full sampling rate (2048 Hz)
        time_full = time  # Already defined earlier using np.linspace(...)

        # Interpolate (resample) the PCA signal at the full resolution of the original signal
        reference_resampled = interp_func(time_full)





        def norm_and_smooth(series, window=6, polyorder=4):
            normed = norm_df(series)
            smoothed = normed
            return pd.Series(smoothed, index=series.index)

        def CalculateCorrCoffs(part,dfdata):
        
            # Normalize the respiratory signal to [0, 1] range
            scaler = MinMaxScaler()
            resp_norm = scaler.fit_transform(part.reshape(-1, 1)).flatten()


            # Plot the normalized respiratory signal
            plt.figure(figsize=(20, 6))
            plt.plot(resp_norm, label="Normalized Respiration", linewidth=4)


            datatime =dfdata['data1'][1:].astype(int)
            
            indexedPCA= resp_norm[datatime]

            plt.plot(datatime,indexedPCA, label="pca at data time")


            dfdataL = norm_and_smooth(dfdata['data2'][1:].astype(float))
            dfdataGCV = norm_and_smooth(dfdata['data3'][1:].astype(float))
            dfdataRGCV = norm_and_smooth(dfdata['data4'][1:].astype(float))
            dfdataCRESO = norm_and_smooth(dfdata['data5'][1:].astype(float))
            dfdataU = norm_and_smooth(dfdata['data6'][1:].astype(float))
            dfdataLFix = norm_and_smooth(dfdata['data7'][1:].astype(float))
            avgResidual = norm_and_smooth(dfdata['data8'][1:].astype(float))
            avgSolution = norm_and_smooth(dfdata['data9'][1:].astype(float))
            # Compute correlation coefficients
            correlations = {
                "L-curve": corr_with_pca(dfdataL, indexedPCA),
                "GCV": corr_with_pca(dfdataGCV, indexedPCA),
                "RGCV": corr_with_pca(dfdataRGCV, indexedPCA),
                "CRESO": corr_with_pca(dfdataCRESO, indexedPCA),
                "U-curve": corr_with_pca(dfdataU, indexedPCA),
                "L-curve fix": corr_with_pca(dfdataLFix, indexedPCA),
                "Avg Residual": corr_with_pca(avgResidual, indexedPCA),
                "Avg Solution": corr_with_pca(avgSolution, indexedPCA)
            }

        # Plot each normalized series
            plt.plot(datatime, dfdataL, label=f"L-curve (r={correlations['L-curve']:.2f})")
            plt.plot(datatime, dfdataGCV, label=f"GCV (r={correlations['GCV']:.2f})")
            plt.plot(datatime, dfdataRGCV, label=f"RGCV (r={correlations['RGCV']:.2f})")
            plt.plot(datatime, dfdataCRESO, label=f"CRESO (r={correlations['CRESO']:.2f})")
            plt.plot(datatime, dfdataU, label=f"U-curve (r={correlations['U-curve']:.2f})")
            plt.plot(datatime, dfdataLFix, label=f"L-curve fix (r={correlations['L-curve fix']:.2f})")
            plt.plot(datatime, avgResidual, label=f"Avg Residual (r={correlations['Avg Residual']:.2f})")
            plt.plot(datatime, avgSolution, label=f"Avg Solution (r={correlations['Avg Solution']:.2f})")

            # Final plot setup

            

            plt.legend()
            plt.xlabel("Time")
            plt.ylabel("Normalized Values")
            plt.title("Comparison of Normalized Regularization Methods")

            plt.grid(True)
            plt.savefig(os.path.join( output_folder, os.path.basename(filepath).split('.')[0]+'.png'), dpi=300, bbox_inches='tight')

            plt.clf()

            return 1

        CalculateCorrCoffs(reference_resampled,dfdata)

    except Exception as e:
        print(f"Skipping {filepath} due to error: {e}")
        continue
