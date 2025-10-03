
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import  lfilter
import pandas as pd
import glob , os

from scipy.io import loadmat ,savemat

# --------- File loading and preprocessing ---------

def find_csv(base_path):
    # Look only in subfolders starting with CHU
    chu_dirs = glob.glob(os.path.join(base_path, "CHU*"))

    matches = []
    for d in chu_dirs:
        # Include only interpSig*.csv
        for p in glob.glob(os.path.join(d, "interpSig*.csv")):
            name = os.path.basename(p)
            # Exclude anything that looks like ECG*.csv just in case
            if not name.lower().startswith("ecg"):
                matches.append(p)

    return matches

base = r"C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical"
paths = find_csv(base)


for filename in paths:
    print('-----------------------------------------------')
    print(f"CSV: {os.path.basename(filename)}")
    try:
        
        # Path to filtered respiration data files and list of bad channels
        # Load the filtered signal from CSV and transpose it
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        df = pd.read_csv(filename, sep=',')
        print(filename)
        

        start = int(1*1000)
   

        Y = df.to_numpy(copy=True)
        Y = Y[ :,start:(-start)]
        print(Y.shape)

        time = np.linspace(0, df.shape[1] / 1000, df.shape[1])
        time = time[start:(-start)]
        print(time.shape)

        # Step 1: Apply FIR filter to each row using lfilter
        # lfilter does not natively broadcast over axis 0, so we use list comprehension or a matrix operation
        h = [0.1, 0.2, 0.0, -0.2, -0.1]
        H = np.array(h)
        Y_filtered = np.apply_along_axis(lambda row: lfilter(H, 1.0, row), axis=1, arr=Y)

        thresholdpercent = 0.5

        windowsize = int(0.50*1000)
        # Step 2: Square the filtered signal
        Y_squared = Y_filtered ** 2
        # Step 3: Moving average (via convolution), applied row-wise
        kernel = np.ones(windowsize) / windowsize
        Y_smoothed = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode='same'), axis=1, arr=Y_squared)
        mean_smoothed = np.mean(Y_smoothed, axis=0)  # Shape: (n_samples,)

        thresholds = np.max(Y_smoothed, axis=1, keepdims=True) * thresholdpercent
        pwm_signals = (Y_smoothed > thresholds).astype(int)

        pwm_signal = (mean_smoothed > np.max(mean_smoothed)*thresholdpercent)


        totalPWM= np.apply_along_axis(lambda coulmn: np.sum(coulmn) > 10,axis=0 ,arr=pwm_signals ).astype(int)

        gradpwm = np.gradient(totalPWM)
        qrsOnset = np.where((gradpwm>0)==True)[0]
        qrsEndset = np.where((gradpwm<0)==True)[0]

        min_gap = 500

        # Keep only the indices that are far enough apart from the previous one
        qrsOnset = qrsOnset[np.insert(np.diff(qrsOnset) > min_gap, 0, True)]
        qrsEndset = qrsEndset[np.insert(np.diff(qrsEndset) > min_gap, 0, True)]

        for _ in range(4):
            if qrsOnset[0] > qrsEndset[0]:  
                qrsEndset = qrsEndset[1:]  # Remove first high peak if it occurs before the first trough
                print("first high_peaks removed")
            # Ensure the signal ends with a high peak
            if qrsOnset[-1] > qrsEndset[-1]:
                print("last low_peaks removed")
                qrsOnset = qrsOnset[:-1]  # Remove last low peak if it occurs after the last high peak
        print('o=o=o=o=o=o=o=o=o=o=o=o=o')
        print(qrsOnset)
        print(qrsEndset)
        print('###############################')
        print(qrsEndset-qrsOnset)
        print('o=o=o=o=o=o=o=o=o=o=o=o=o')
        
        qrsRange = np.mean(qrsEndset-qrsOnset)+100
        beatRange = np.mean(qrsOnset[1:]-qrsOnset[:-1])+100



        plt.figure(figsize=(12, 5))

        # Plot all original signals (absolute value) in gray
        for i in range(Y.shape[0]):
            plt.plot(time, np.abs(Y[i,:]), alpha=0.2, color='gray')

        # Plot mean smoothed
        plt.plot(time, (mean_smoothed*(Y.max())/(mean_smoothed.max()*2)) , label='Mean Smoothed scaled', color='blue')

        # Plot PWM signal
        plt.vlines(time[qrsOnset],ymin=plt.ylim()[0], ymax=plt.ylim()[1], label='Qrsonset', color='green')

        # Plot gradient of PWM
        plt.vlines(time[qrsEndset],ymin=plt.ylim()[0], ymax=plt.ylim()[1],  label='Qrsendset', color='red')

        plt.title(f"{filename}")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.grid(True)
        plt.legend(loc=2)
        plt.tight_layout()
        plt.savefig(f"C:\\Users\\emir.ege-nemutlu\\Desktop\\resp\\Clinical\\qrspng\\{(os.path.basename(filename)).removesuffix(".csv")}_plot.png")
        plt.close()


        qrsOnset = qrsOnset + start
        struct_fields = {
            'Beat_range': beatRange.astype('double'),
            'QRS_range': qrsRange.astype('double'),

            'QRS_Onset_auto': (qrsOnset).reshape(-1, 1).astype('double') 
        }

        # Save to the new .mat file


                
     
        base_folder = os.path.dirname(filename)  # Go up from Signals_Beats
        save_dir = os.path.join(base_folder, "Signal_select_auto")
        # os.makedirs(save_dir, exist_ok=True)
        new_matfile_path = os.path.join(save_dir, (os.path.basename(filename).removesuffix(".csv")+'.mat'))
        os.makedirs(os.path.dirname(new_matfile_path), exist_ok=True)
        savemat(new_matfile_path, {'Signal_select_Beat': struct_fields})

        print(len(qrsOnset[:-1]))
        print(len(qrsEndset[1:]))
        print(np.mean(qrsEndset-qrsOnset))

        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
    
    except Exception as e:
        print(f"Skipping {filename} due to error: {e}")
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        continue


