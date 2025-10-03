clc
clear all
close all

dossier = "U:\Programmes\CardioResp\Pig1\EpiPacingSite-223.mat";  % Path to the ECG .mat file
dossier2 = "U:\Programmes\CardioResp\Pig1\Filtered\";              % Output directory for filtered signals
name = "EpiPacingSite-223";                                       % Base name for saved output
fstruct = dir("U:\Programmes\CardioResp\Pig1\EpiPacingSite-133.mat");  % File list 

for m = 2 % Loop index 
    file = dossier;                      % Get file path
    load(file)                           % Load the MAT file containing 'rec' and possibly 'info'

    if(m > length(fstruct))             % This block would only execute if m exceeds number of files
        ECG = signalObject([rec.vest.Ve],'samp',2048);           % Create ECG object with sampling rate
        ECG.applyLowPassFilter(60,10,'FIR');                     % Apply low-pass filter at 60 Hz, order 10
        temp = ECG.processedSignal;                              % Store result in temp
    else
        ECG = signalObject([rec.vest.Ve],'samp',2048);           % Create ECG object from raw signal
        ECG.applyLowPassFilter(60,6,'FIR');                      % Apply 60 Hz low-pass filter with order 6
        ECG.removeBaseline('savitzkygolay', [3 5000]);           % Remove baseline using Savitzky-Golay filter
        sig = ECG.processedSignal;                               % Extract processed ECG signal
        [locs,lengthSpike] = findSpikes(sig,2);                  % Detect pacing spikes and their duration

        ECG = signalObject([rec.vest.Ve_filtered],'samp',2048);  % Reload filtered version if available
        ECG.applyLowPassFilter(60,6,'FIR');                      % Apply low-pass filter again
        temp = ECG.processedSignal;                              % Store the signal to be cleaned

        for n = 1:length(locs)                                   % Loop through each detected spike
            x = [locs(n)-70:locs(n)-round(lengthSpike/2)-1, locs(n)+round(lengthSpike/2)+1:locs(n)+70];  % Exclude spike center

            if(find(x<1))                                        % Remove indices below start of signal
                x(find(x<1))=[];
            end
            if(find(x>size(temp,2)))                             % Remove indices beyond signal length
                x(find(x>size(temp,2)))=[];
            end

            y = ECG.processedSignal(:,x);                        % Extract signal around spike (excluding center)
            xx = x(1):x(end);                                    % Create full range for interpolation
            yy = spline(x,y,xx);                                 % Interpolate using spline to fill the spike
            temp(:,xx)=yy;                                       % Replace spike region with interpolated values
        end

        figure(2); clf
        plot(rec.vest.Ve_filtered(14,:))                         % Plot original filtered signal (channel 14)
        hold on
        plot(temp(14,:))                                         % Overlay with spike-cleaned version
    end

    ECG = signalObject([temp],'samp',2048);                      % Create ECG object from cleaned signal
    ECG.applyNotchFilter(50,2,'BandStop');                       % Remove 50 Hz powerline noise
    ECG.applyNotchFilter(100,2,'BandStop');                      % Remove 100 Hz harmonic
    ECG.applyNotchFilter(150,2,'BandStop');                      % Remove 150 Hz harmonic
    ECG.removeBaseline('savitzkygolay', [3 5000]);               % Final baseline correction

    figure(4); clf
    plot(ECG.rawSignal(14,4000:6000),'k')                        % Plot raw signal (black)
    hold on
    plot(rec.vest.Ve(14,4000:6000),'b')                          % Original unprocessed signal (blue)
    plot(ECG.processedSignal(14,4000:6000),'r')                  % Final filtered and cleaned signal (red)

    oldFilter = rec.vest.Ve;                                     % Store original signal
    rec.vest.Ve_filtered = ECG.processedSignal;                  % Update struct with filtered version

    ECG = signalObject([rec.vest.Ve],'samp',2048);               % Reload original signal for comparison
    ECG.applyNotchFilter(50,2,'BandStop');
    ECG.applyNotchFilter(100,2,'BandStop');
    ECG.applyNotchFilter(150,2,'BandStop');
    ECG.removeBaseline('savitzkygolay', [3 5000]);

    figure(3); clf
    subplot(2,2,1)
    plot(rec.vest.Ve(:,4000:6000)')                              % Raw signal
    ylim([-3 3])
    subplot(2,2,2)
    plot(oldFilter(:,4000:6000)')                                % Previously saved filtered version
    ylim([-3 3])
    subplot(2,2,3)
    plot(ECG.processedSignal(:,4000:6000)')                      % Cleaned after full processing
    ylim([-3 3])
    subplot(2,2,4)
    plot(rec.vest.Ve_filtered(:,4000:6000)')                     % Final version saved in struct
    ylim([-3 3])

    save([dossier2 + name + "_filtered"],'rec','info')           % Save final struct
    csvwrite(dossier2+ name + "_filtered.csv",rec.vest.Ve_filtered)  % Save filtered signal as CSV
end
