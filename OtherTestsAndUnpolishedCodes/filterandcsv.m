dossier2 = "C:\Users\emir.ege-nemutlu\Desktop\resp"; 
name = "EndoPacingSite-LV3";
writematrix(rec.vest.Ve_filtered, fullfile(dossier2, name + "_filtered.csv"));  % Save filtered signal as CSV

data = rec.vest.Ve_filtered;

% Sampling frequency (you should set this properly)
Fs = 2048;  % Example: 1000 Hz
t = (0:size(data, 2)-1) / Fs;  % Time vector in seconds

