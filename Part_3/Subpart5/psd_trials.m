function [standard_psd, psd50, psd150, samples_50, samples_150] = psd_trials(RRI_trial, fs)
    addpath(genpath([fileparts(pwd), '\Subpart1']));
    
    trial = zscore(RRI_trial);
    
    samples_50 = fs*50;
    samples_150 = fs*150;
    
    n_segments50 = floor(length(trial)/samples_50);
    n_segments150 = floor(length(trial)/samples_150);
    
    % standard PSD of the trial
    standard_psd = pgm(trial);
    
    % averaging for window size of 50
    psd50 = zeros(1, samples_50);
    
    for i = 0:n_segments50-1
        window = trial((samples_50*i)+1 : samples_50*(i+1));
        psd50 = psd50 + (pgm(window) ./ n_segments50);
    end

    % averaging for window size of 150
    psd150 = zeros(1, samples_150);
    
    for i = 0:n_segments150-1
        window = trial((samples_150*i)+1 : samples_150*(i+1));
        psd150 = psd150 + (pgm(window) ./ n_segments150);
    end
end