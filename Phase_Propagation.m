%% Phase Propagation

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
    'dom_freq', {}, 'mean_psd', {}, 'filt_mean_psd', {}, 'LFP_bandpass', {});

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

    cleanData(i).LFP_bandpass = zeros(size(cleanData(i).lfp));

end

%%

domFreq = 12.5;
[b, a] = butter(2, [domFreq-1 domFreq+1]/(fs/2), "bandpass");

figure
freqz(b, a, [], fs)
xlim([0 30])
subplot(2,1,2)
xlim([0 30])
set(gcf, 'WindowState', 'maximized')
if saveFig
    saveas(gcf,'pics_p2/01_filter.png');
end

for i = 1:numCh
    [b, a] = butter(2, [cleanData(i).dom_freq-1 cleanData(i).dom_freq+1]/(fs/2), "bandpass");
    for j = 1:numTrials
    cleanData(i).LFP_bandpass(:,j) = filtfilt(b, a, cleanData(i).lfp(:,j));
    end
end

% plotting a sample trial and it's filtered
figure
plot(tVec, cleanData(9).lfp(:,67), 'Color',	"#77AC30", 'LineWidth', 1.5)
hold on
plot(tVec, cleanData(9).LFP_bandpass(:,67) , 'Color', "#D95319", 'LineWidth', 1.5)
grid on
ylabel('Amplitude [\muV]')
xlabel('Time [s]')
legend('Original','Filtered')
title('LFP before and after Filtering')
set(gcf, 'WindowState', 'maximized')
if saveFig
    saveas(gcf,'pics_p2/02_filtered.png');
end

%% PART B
phi = nan*zeros(numCh, numTrials, length(tVec));

for i = 1:numCh
    for j = 1:numTrials
        phi(i, j, :) = angle(hilbert(cleanData(i).LFP_bandpass(:,j)));
    end
end

% plotting phase of a sample trial
figure
plot(tVec, reshape(phi(6,89,:), 1, 641) , 'Color', "#D95319", 'LineWidth', 1.5)
grid minor
ylabel('Phase [radian]')
xlabel('Time [s]')
title('Phase of LFP Signal')
set(gcf, 'WindowState', 'maximized')
if saveFig
    saveas(gcf,'pics_p2/03_phase.png');
end

%% PART C
trl = 123;

cosPhiMat1 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    cosPhiMat1(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = cos(phi(i, trl, :));
end

figure
tiledlayout(4,4);
for i = 500:500+15
    nexttile
    imagesc(1:10 ,0:5 ,cosPhiMat1(:,:,i));
    colormap hot
    clim([-1 1])
    title(['time = ', num2str(tVec(i)), 's'])
end
c = colorbar('location', 'manual');
set(c, 'position', [0.95 0.1 0.02 0.8]);
set(gcf, 'WindowState', 'maximized')
sgtitle(['Wave Prop Snapshot of Trial ', num2str(trl)])

if saveFig
    saveas(gcf,'pics_p2/04_snapshot.png');
end

% extracting phase of 3 other trials
cosPhiMat2 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    cosPhiMat2(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = cos(phi(i, 66, :));
end

cosPhiMat3 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    cosPhiMat3(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = cos(phi(i, 400, :));
end

cosPhiMat4 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    cosPhiMat4(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = cos(phi(i, 263, :));
end


% making video demo
frames = [];
figure
set(gcf, 'WindowState', 'maximized')
for i = 1:length(tVec)
    subplot(2,2,1)
    imagesc(1:10 ,0:5 ,cosPhiMat1(:,:,i));
    colormap hot
    clim([-1 1])
    colorbar
    title(['Demo for Trial', num2str(trl), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    subplot(2,2,2)
    imagesc(1:10 ,0:5 ,cosPhiMat2(:,:,i));
    colormap hot
    clim([-1 1])
    colorbar
    title(['Demo for Trial', num2str(66), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    subplot(2,2,3)
    imagesc(1:10 ,0:5 ,cosPhiMat3(:,:,i));
    colormap hot
    clim([-1 1])
    colorbar
    title(['Demo for Trial', num2str(400), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    subplot(2,2,4)
    imagesc(1:10 ,0:5 ,cosPhiMat4(:,:,i));
    colormap hot
    clim([-1 1])
    colorbar
    title(['Demo for Trial', num2str(263), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    frames = [frames getframe(gcf)];
end

writerObj = VideoWriter("Wave_Propagation");
writerObj.FrameRate = fs/10;
writerObj.Quality = 100;
open(writerObj);
for i=1:length(frames)
    frame = frames(i) ;
    writeVideo(writerObj,frame);
end
close(writerObj)

%% PART D

PhiMat1 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    PhiMat1(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = phi(i, trl, :);
end

[gx, gy] = gradient(unwrap(PhiMat1));
propDir1 = cell(1, length(tVec));

for i = 1:length(tVec)
    propDir1{i} = atan2(gy(:,:,i), gx(:,:,i));
end

figure
polarhistogram(cell2mat(propDir1), 60, 'FaceColor', "#77AC30");
title(['Propagation Direction Histogram of Trial ', num2str(trl)])
set(gcf, 'WindowState', 'maximized')
if saveFig
    saveas(gcf,'pics_p2/05_prop_dir_all.png');
end

figure
tiledlayout(4,4);
for i = 500:500+15
    nexttile
    polarhistogram(propDir1{i}, 30, 'FaceColor', "#77AC30");
    title(['Time=', num2str(tVec(i)), ' s'])
end
set(gcf, 'WindowState', 'maximized')
sgtitle(['Propagation Direction Histogram of Trial ', num2str(trl), ' at some samples'])
if saveFig
    saveas(gcf,'pics_p2/06_prop_dir_snapshot.png');
end

% extracting prop direction of 3 other trials
PhiMat2 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    PhiMat2(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = phi(i, 66, :);
end
[gx, gy] = gradient(unwrap(PhiMat2));
propDir2 = cell(1, length(tVec));
for i = 1:length(tVec)
    propDir2{i} = atan2(gy(:,:,i), gx(:,:,i));
end

PhiMat3 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    PhiMat3(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = phi(i, 400, :);
end
[gx, gy] = gradient(unwrap(PhiMat3));
propDir3 = cell(1, length(tVec));
for i = 1:length(tVec)
    propDir3{i} = atan2(gy(:,:,i), gx(:,:,i));
end

PhiMat4 = nan*zeros(5, 10, length(tVec));
for i = 1:numCh
    PhiMat4(cleanData(i).ch_position_y, cleanData(i).ch_position_x, :) = phi(i, 263, :);
end
[gx, gy] = gradient(unwrap(PhiMat4));
propDir4 = cell(1, length(tVec));
for i = 1:length(tVec)
    propDir4{i} = atan2(gy(:,:,i), gx(:,:,i));
end

%% PART E
% making video demo
frames = [];
figure
set(gcf, 'WindowState', 'maximized')
for i = 1:length(tVec)
    subplot(2,2,1)
    polarhistogram(propDir1{i}, 30, 'FaceColor', "#77AC30");
    title(['Demo for Trial ', num2str(123), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    subplot(2,2,2)
    polarhistogram(propDir2{i}, 30, 'FaceColor', "#77AC30");
    title(['Demo for Trial ', num2str(66), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    subplot(2,2,3)
    polarhistogram(propDir3{i}, 30, 'FaceColor', "#77AC30");
    title(['Demo for Trial ', num2str(400), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    subplot(2,2,4)
    polarhistogram(propDir4{i}, 30, 'FaceColor', "#77AC30");
    title(['Demo for Trial ', num2str(263), ' at Time = ', num2str(tVec(i), '%.2f'), ' s'])

    frames = [frames getframe(gcf)];
end

writerObj = VideoWriter("Wave_Propagation_Direction");
writerObj.FrameRate = fs/10;
writerObj.Quality = 100;
open(writerObj);
for i=1:length(frames)
    frame = frames(i) ;
    writeVideo(writerObj,frame);
end
close(writerObj)

%% PART F

%% PART G
domFreq = 12.5;
omega = 2*pi*domFreq;
tmp = zeros(length(tVec), 5, 10);
for i = 1:length(tVec)
    for j = 1:5
        tmp(i, j, :) = reshape(PhiMat1(j, :, i), 1, 10);
    end
end
deltaT = gradient(tmp) / omega;

distance_elec = 4e-04;

v_wave = distance_elec./deltaT;

figure
edges = linspace(0, 2, 50);
histogram(abs(v_wave), edges);
xlim([0 2])
grid minor
xlabel('Wave Speeds [m/s]')
ylabel('Count')
title('Wave Propagating Speed Histogram')
set(gcf, 'WindowState', 'maximized')
if saveFig
    saveas(gcf,'pics_p2/07_wave_speed.png');
end















