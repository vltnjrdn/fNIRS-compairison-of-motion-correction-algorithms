%%%% COMPAIRISON OF MOTION CORRECTION ALGORITHMS %%%%
%%%%%%%%%% ONE TIME WINDOW AND ONE CHANNEL %%%%%%%%%%

%% README
% I kept this skript because some of the algorithms show different
% performace depending on wheter they are applied to the whole dataset
% compared to a specific time window

% for this skript to run as it is right now you will need:
% 1. homer2
% 2. one fNIRS dataset from the taste experiment (I used 2025-11-13_001)
% 3. the compute_PSD_oxyHb.m function I wrote to get a table with the Power 
%    of the signal in dB in 4 different frequency bands 
%    ('0.015–0.15 Hz','0.2–0.4 Hz','0.7–1.3 Hz','2–5 Hz')

%% motion correction algorithms used
    % 1. spline interpolation (hmrMotionCorrectSpline) + SG (hmrMotionCorrectSG)
    % 2. spline interpolation (hmrMotionCorrectSpline) + RLOESS (hmrMotionCorrectRLOESS)
    % 3. spline interpolation (hmrMotionCorrectSpline) + Wavelet Correction (hmrMotionCorrectWavelet)
    % 4. spline interpolation (hmrMotionCorrectSpline) + PCA (hmrMotionCorrectPCA_Ch)
    % (Jahani et al., 2018)

%% evaluation metrics
    % 1. Number of motion artifacts (n_motion_artifacts)
        % The Number of detected motion artifacts is expected to lower
        % after motion correction. For detecting motion artifacts we'll
        % apply the hmrMotionArtifact and/or the hmrMotionArtifactByChannel
        % function to the OD time Series (Pinto et al., 2023 and Brigadoi et al., 2014)
        % we will count the timepoints labelled as motion artifacts in the
        % whole dataset and only those parts that were labelled as motion
        % artifacts in the uncorrected dataset.
    % 2. Motion ratio (motion_ratio)
        % ratio between the cummulative time of segments in the data that
        % are considered to be motion artifacts by the hmrMotionArtifactByChannel 
        % function to the total acquisition time. This will also be
        % calculated in the whole dataset but also only the parts labelled
        % as motion artifacts before.
        % (von Lühmann et al., 2020)
    % 3. Power Spektral Density (psd)
        % We expect the Power (dB) of physiological frequencies 
        % (respiration (0.2-0.4 Hz), heartbeat (0.7-1.3 Hz)) and cerebral 
        % signal (0.015-0.15 Hz) to get stronger and the non-physiological 
        % associated with noise (2-5 Hz) to get weaker. This measure is 
        % taken before applying the bandpass filter.
    % 4. mean amplitude: (mean_A)
        % Generally, the mean amplitude is expected to be lower after
        % motion correction but we are more interested in compairing the 
        % mean amplitude between the different motion correction algorithms
        % to check for oversmoothing. The mean amplitude of each signal was
        % calculated as the root-mean-square (RMS) value, defined as the 
        % square root of the mean of the squared signal samples
        % (Gemignani and Gervain, 2021)
    % 5. Visual inspection of the plotted HbO and HbR signal

%% workflow
% 1. Import data
% 2. Convert Intensity into OD (hmrIntensity2OD) and extract experiment data
% 3. Motion artifact correction
    % 0 = no correction (noCorr)
    % 1 = pipeline 1 (SG)
    % 2 = pipeline 2 (RLOESS)
    % 3 = pipeline 3 (Wavelett)
    % 4 = Pipeline 4 (PCA)
% 4. Motion artifact detection (n_motion_artifacts, motion_ratio, decline)
% 5. Convert OD into Concentration (hmrOD2Conc)
% 6. Compute PSD
% 7. Bandpass filter (hmrBandpassFilt, 0.01, 0.5)
% 8. Plotting of HbO and HbR
% 9. Mean Amplitude (sqrt(mean(data(:).^2)))

%% 1. Import data
data = import_homer2('2025-11-13_001.nirs');

% parameters
SamplingRate = 1/t(2,1); % Sampling Rate
duration = 28; % duration of the experiment in minutes

stim1 = find(s(:,1)==1); % taste delivery
stim2 = find(s(:,2)==1); % hold taste
stim3 = find(s(:,3)==1); % swallow and rinse
stim4 = find(s(:,4)==1); % reinsert mouthpiece & rest
stim5 = find(s(:,5)==1); % continuing shortly
expStart = stim1(1); % start of the experiment
expEnd = stim1(1) + duration*60*SamplingRate; % end of the experiment

%% User-defined parameters
% choose the channel
channel_to_use = 1;
% choose the time indices
t_start = stim1(1); 
t_end = stim1(2);

time = (t_start:t_end)/SamplingRate; % user defined time vector
time_plot = time/60; % time in minutes for plotting

algNames = {'noCorr','SG','RLOESS','Wavelett','PCA'}; % names of the motion correction algorithms
nDatasets = length(algNames); % number of datasets

%% 2. Convert intensity to optical density and extract time window
dod = hmrIntensity2OD(d);

% extract time window
dod = dod(t_start:t_end,:);

%% 3. motion artifact correction
% (takes abt 40s, depends on time window size)

% tIchCh for spline interpolation and PCA that only correct segments marked as motion artifact
tIncMan = ones(size(d,1),1);
SD.MeasListAct = true(size(SD.MeasList,1),1);
[tInc, tIncCh] = hmrMotionArtifactByChannel(dod, time, SD, tIncMan, 0.5, 1, 15, 0.2);

% 0 = no correction (noCorr)
dod_0 = dod;

% spline interpolation
dod_spline = hmrMotionCorrectSpline(dod, time, SD, tIncCh, 0.99);

% 1 = savitzky golay (SG)
dod_1 = hmrMotionCorrectSG(dod_spline, time, 10, 1);

% 2 = rloess (RLOESS)
dod_2 = hmrMotionCorrectRLOESS(dod_spline, time, 0.02, 1);

% 3 = Wavelet correction (Wavelett)
dod_3 = hmrMotionCorrectWavelet(dod_spline, SD, 1.5);

% 4 = PCA
dod_4 = hmrMotionCorrectPCA_Ch( SD, dod_spline, tInc, tIncCh, 0.9);

% save corrected datasets in a list
dod_list = {dod_0, dod_1, dod_2, dod_3, dod_4};

%% 4. Motion artifact detection

% prepare table for measures
rowNames = {'motion artifacts (time window)', 'motion ratio (time window)', 'motion artifacts (motion mask)', 'motion ratio (motion mask)'};
eval_motion_artifacts = array2table(NaN(length(rowNames), length(algNames)), 'VariableNames', algNames, 'RowNames', rowNames);

% compute original motion artifact mask before correction
[~, tIncCh_orig] = hmrMotionArtifactByChannel(dod_0, time, SD, tIncMan, 0.5, 1, 15, 0.2);
bad_mask = (tIncCh_orig == 0); % mask for parts marked as motion artifact (flag == 0)

% loop over each corrected dataset
for idx = 1:length(dod_list)

    dod = dod_list{idx};
    
    % motion artifact masks after correction
    [tInc, tIncCh] = hmrMotionArtifactByChannel(dod, time, SD, tIncMan, 0.5, 1, 15, 0.2);

    % number and ratio of motion artifacts after correction for whole dataset
    n_all = sum(tIncCh(:,channel_to_use) == 0); % number
    r_all = n_all / (numel(tIncCh(:,channel_to_use))); % ratio

    eval_motion_artifacts{'motion artifacts (time window)', algNames{idx}} = n_all;
    eval_motion_artifacts{'motion ratio (time window)', algNames{idx}} = r_all;

    % number and ratio of motion artifacts after correction only inside the
    % parts already marked as motion artifacts before
    
    n_motion = sum(tIncCh(bad_mask(:,channel_to_use)) == 0); % number 
    r_motion = n_motion / numel(bad_mask(:,channel_to_use)); % ratio 

    eval_motion_artifacts{'motion artifacts (motion mask)', algNames{idx}} = n_motion;
    eval_motion_artifacts{'motion ratio (motion mask)', algNames{idx}} = r_motion;

end

disp(eval_motion_artifacts)

%% 5. Convert OD into Concentration and extract oxygenated and deoxygenated Haemoglobin

% prepare list for the concentrations
dc_list = cell(size(dod_list));

for idx = 1:length(dod_list)
    dc_list{idx} = hmrOD2Conc(dod_list{idx}, SD, [6 6]); 
    % dc_list has 3 dimensions: time x chromophores(3) x channels(28) 
end

% Extract oxygenated Haemoglobin
HbO_list = cell(size(dc_list));

for idx = 1:length(dc_list)
    HbO_list{idx} = squeeze(dc_list{idx}(:, 1, channel_to_use));  % time x chromophores(3) x channel(1)
end

% Extract desoxygenated Haemoglobin
HbR_list = cell(size(dc_list));
for idx = 1:length(dc_list)
    HbR_list{idx} = squeeze(dc_list{idx}(:, 2, channel_to_use));  % time x chromophores(3) x channel(1)
end


%% 6. Compute PSD

% compute PSD for four different frequency bands (in dB) with compute_PSD_oxyH function
eval_PSD = compute_PSD_oxyHb(HbO_list, SamplingRate, algNames);
disp(eval_PSD);


% find min and max of PSD and frequency to adjust xlim and ylim
allPSD = [];
allF   = [];

for col = 1:nDatasets
    data = HbO_list{col};
    [Pxx, f] = pwelch(data, [], [], [], SamplingRate);
    allPSD = [allPSD; 10*log10(Pxx(:))];
    allF   = f;
end

yLim = [min(allPSD) max(allPSD)];
xLim = [min(allF) max(allF)];

% plot the channel
figure('Name','PSD Comparison - Selected Channel');

for col = 1:nDatasets
    data = HbO_list{col};
    [Pxx, f] = pwelch(data, [], [], [], SamplingRate);

    subplot(1, nDatasets, col);
    plot(f, 10*log10(Pxx));
    grid on;

    xlim(xLim);
    ylim(yLim);

    title(['Ch ', num2str(channel_to_use), ' - ', algNames{col}]);

    if col == 1
        ylabel('Power/Frequency (dB)');
    end
    xlabel('Frequency (Hz)');
end

%% 7. Bandpass filter

% prepare list for filtered signal
HbO_filt_list = cell(size(HbO_list));
HbR_filt_list = cell(size(HbR_list));

% apply bandpass filter to oxygenated and deoxygenated signal
for idx = 1:length(HbO_list)
    HbO_filt_list{idx} = hmrBandpassFilt(HbO_list{idx}, time, 0.01, 0.5);
    HbR_filt_list{idx} = hmrBandpassFilt(HbR_list{idx}, time, 0.01, 0.5);
end

% if you don't want to apply any bandpass filter:
% HbO_filt_list = HbO_list;
% HbR_filt_list = HbR_list;

%% 8.2 plot a specific time sequence

% plot all algorithms in one plot
figure;
hold on;
for idx = 1:length(HbO_filt_list)
    data = HbO_filt_list{idx};
    plot(time_plot, data, 'LineWidth', 1.5);

    xline((time_plot(1) + 0.1/60),  '--k', 'Taste delivery', 'LabelVerticalAlignment', 'top', 'FontSize', 10, 'HandleVisibility','off');
    xline((time_plot(1) + 21/60), '--k', 'Swallow & rinse', 'LabelVerticalAlignment', 'top', 'FontSize', 10, 'HandleVisibility','off');
    xline((time_plot(1) + 41/60), '--k', 'Rest', 'LabelVerticalAlignment', 'top', 'FontSize', 10, 'HandleVisibility','off');
end

hold off;
grid on;

title(['HbO after different corrections - Ch ', num2str(channel_to_use)]);
xlabel('Time (min)');
ylabel('\DeltaHbO (\muM)');
xlim([time_plot(1) time_plot(end)]);

legend(algNames, 'Location','best');

% plot HbO and HbR for all algorithms in different plots

% global min and global max
allData = [];
for idx = 1:nDatasets
    allData = [allData; HbO_filt_list{idx}; HbR_filt_list{idx}];
end
yMin = min(allData);
yMax = max(allData);

% HbO and HbR
colors = lines(2);
figure('Name','HbO/HbR for every algorithm');

for idx = 1:nDatasets
    subplot(2,3,idx);
    hold on;
    
    % HbO
    dataO = HbO_filt_list{idx};
    plot(time_plot, dataO, 'Color', colors(2,:), 'LineWidth', 1.5);

    % HbO
    dataR = HbR_filt_list{idx};
    plot(time_plot, dataR, 'Color', colors(1,:), 'LineWidth', 1.5);
    
    xline((time_plot(1) + 0.1/60),  '--k', 'Taste delivery', 'LabelVerticalAlignment', 'top', 'FontSize', 8, 'HandleVisibility','off');
    xline((time_plot(1) + 21/60), '--k', 'Swallow & rinse', 'LabelVerticalAlignment', 'top', 'FontSize', 8, 'HandleVisibility','off');
    xline((time_plot(1) + 41/60), '--k', 'Rest', 'LabelVerticalAlignment', 'top', 'FontSize', 8, 'HandleVisibility','off');
        
    hold off;
    grid on;

    title(['Ch ', num2str(channel_to_use), ' - ', algNames{idx}]);
    xlabel('Time (min)');
    ylabel('\DeltaHb (\muM)');
    
    ylim([yMin yMax]);
    xlim([time_plot(1) time_plot(end)]);

    legend({'\DeltaHbO', '\DeltaHbR'}, 'Location','best');
end

%% mean amplitude

% Prepare table: 1 row per dataset
eval_meanAmp = table('Size',[1, nDatasets], ...
                      'VariableTypes', repmat("double",1,nDatasets));

% Fill table
for idx = 1:nDatasets
    data = HbO_filt_list{idx};      % time
    eval_meanAmp{1,idx} = sqrt(mean(data(:).^2));  % mean over all time points and channels
end

% Set column names
eval_meanAmp.Properties.VariableNames = algNames;

% Display
disp(eval_meanAmp)
