function PSD_table = compute_PSD_oxyHb_whole(oxyHb_list, fs, algNames)
% Computes PSD metrics for oxygenated hemoglobin

% Frequency bands
freq_bands = [0.015 0.15; % brain-related HbO
              0.2  0.4;   % respiration
              0.7  1.3;   % cardiac signal
              2    5];    % high-frequency noise

band_names = {'0.015–0.15 Hz (cerebral signal)','0.2–0.4 Hz (respiration)','0.7–1.3 Hz (heartbeat)','2–5 Hz (noise)'};

% number of datasets
nDatasets = length(oxyHb_list);

% prepare the table
PSD_table = table('Size',[length(band_names) nDatasets], ...
                  'VariableTypes', repmat("double",1,nDatasets), ...
                  'RowNames', band_names, ...
                  'VariableNames', algNames);

% loop over datasets and compute PSD
for idx = 1:nDatasets
    data = oxyHb_list{idx};  % time × channels
    nCh = size(data,2);
    
    psd_ch = zeros(nCh, size(freq_bands,1));
    
    for ch = 1:nCh
        x = data(:,ch);
        
        if any(~isnan(x))
            [Pxx, F] = pwelch(x, [], [], [], fs);
            for b = 1:size(freq_bands,1)
                idx_band = F >= freq_bands(b,1) & F <= freq_bands(b,2);
                psd_ch(ch,b) = mean(Pxx(idx_band));
            end
        else
            psd_ch(ch,:) = NaN;
        end
    end
    
    % Average across channels and convert to dB
    PSD_table{:,idx} = 10*log10(mean(psd_ch,1,'omitnan')' + eps);
end
end
