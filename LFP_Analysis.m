%% LFP Analysis

% Amirreza Bahramani
% Advanced Neuroscience
% April 2023
% HW 4
% Travelling Wave

close all;
clear;
clc;

%% LOADING
load('ArrayData.mat')
load('CleanTrials.mat')

%% PARAMETERS
saveFig = true;

numCh = 48;
numTrials = 490;
fs = 200;
pwelchWinSize = 256;
cleanData = struct('ch_position_y', {}, 'ch_position_x', {}, 'lfp', {}, ...
    'dom_freq', {}, 'mean_psd', {}, 'filt_mean_psd', {});

%% PART A
% finding the most dominant frequency 

tVec = Time;

for i = 1:numCh
    [cleanData(i).ch_position_y, cleanData(i).ch_position_x] = find(ismember(ChannelPosition,i));

    tmp = chan(i).lfp;
    cleanData(i).lfp = tmp(:, Intersect_Clean_Trials);

    tmp = reshape(cleanData(i).lfp, 1, []);
    [cleanData(i).mean_psd, fVec] = pwelch(tmp, pwelchWinSize, [], [], fs);
    pinkNoise = 1./fVec;
    pinkNoise(1) = pinkNoise(2);
    cleanData(i).mean_psd = 10*log10(cleanData(i).mean_psd);
%%%%%%%%%%%%%%%%%%%
    [~,~,cleanData(i).filt_mean_psd] = regress(cleanData(i).mean_psd, pinkNoise);
%%%%%%%%%%%%%%%%%%%
    [~, idx] = max(cleanData(i).filt_mean_psd); % index of maximum PSD value
    cleanData(i).dom_freq = fVec(idx);

    cleanData(i).wavelet = zeros(55, 641);

end

% Plotting figures
figure;
hold on
for i = 1:numCh
    plot(fVec, cleanData(i).mean_psd)
end
ylabel('PSD [dB/Hz]')
xlabel('Frequency [Hz]')
grid on
title('Mean Power Spectrum of All Trials in each Channel before Denoising')
if true
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,'pics/01_PSD_before_denoise.png');
end

figure;
hold on
for i = 1:numCh
    plot(fVec, cleanData(i).filt_mean_psd)
end
ylabel('PSD [dB/Hz]')
xlabel('Frequency [Hz]')
grid on
title('Mean Power Spectrum of All Trials in each Channel after Denoising')
if saveFig
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,'pics/02_PSD_after_denoise.png');
end

%% PART B
domFreqsMat = nan*zeros(5,10);
for i = 1:numCh
    domFreqsMat(cleanData(i).ch_position_y, cleanData(i).ch_position_x) = cleanData(i).dom_freq;
end

figure
h = heatmap(1:10, 1:5, domFreqsMat);
h.Title = 'Dominant Frequency of each Electrode [Hz]';
h.XLabel = 'Horizontal Position';
h.YLabel = 'Vertical Position';
set(gcf, 'WindowState', 'maximized')
if true
    saveas(gcf,'pics/03_elec_matrix.png');
end

%% PART C
figure
tiledlayout(5,10);
for i = 1:numCh
    for j = 1:numTrials
        [cfs,frq] = cwt(cleanData(i).lfp(:,j), fs, FrequencyLimits=[1 45]);
        cleanData(i).wavelet = cleanData(i).wavelet + abs(cfs);
    end
    cleanData(i).wavelet = cleanData(i).wavelet/numTrials;
    
    nexttile
    surface(tVec,frq,abs(cleanData(i).wavelet))
    axis tight
    shading flat
    xlabel("Time [s]")
    ylabel("Frequency [Hz]")
    colorbar
    % xline(0, 'LineWidth', 2, 'Color', "#A2142F")
    title("LFP Spectrogram of Channel ", num2str(i))
end
set(gcf, 'WindowState', 'maximized')
if true
    saveas(gcf,'pics/04_all_spectrogram.png');
end


mean_wt = zeros(55,641);
for i = 1:numCh
    mean_wt = mean_wt + abs(cleanData(i).wavelet);
end
mean_wt = mean_wt/numCh;

figure
surface(tVec,frq,abs(mean_wt))
axis tight
shading flat
xlabel("Time [s]")
ylabel("Frequency [Hz]")
colorbar
xline(0, 'LineWidth', 2, 'Color', "#A2142F")
title("LFP Spectrogram of All Channels")
set(gcf, 'WindowState', 'maximized')
if true
    saveas(gcf,'pics/05_mean_spectrogram.png');
end

%% PART D

% In report