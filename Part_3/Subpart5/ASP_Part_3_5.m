clear all
close all
clc

%% Set up
pwdd = erase(fileparts(pwd), 'Part_3');
addpath(genpath([pwdd, 'Part_2\Subpart5']));

load RRI_2_5.mat

which = [1 2 3];

full_signal = [RRI_trial1 RRI_trial2 RRI_trial3];
full_signal = zscore(full_signal);
plot(full_signal)

trials = {RRI_trial1; RRI_trial2; RRI_trial3};
frequencies = [RRI_fs1 RRI_fs2 RRI_fs3];
clear RRI_trial1 RRI_trial2 RRI_trial3 RRI_fs1 RRI_fs2 RRI_fs3

if (which(1) >= 1 && which(1) <=3)
    for i = which
        [psd, psd50, psd150, n50, n150] = psd_trials(trials{i}, frequencies(i));

        figure()
        plot(0:1/length(psd):1-1/length(psd), psd)
        title(sprintf('Standard periodogram for Trial %.0f', i))
        xlabel('normalised frequency (2\pi rad/sample)')
        ylabel('Magnitude (Au)')
        xlim([0 0.2])

        figure()
        plot(0:1/n50:1-1/n50, psd50)
        title(sprintf('Averaged periodogram for 50s of Trial %.0f', i))
        xlabel('normalised frequency (2\pi rad/sample)')
        ylabel('Magnitude (Au)')
        xlim([0 0.2])

        figure()
        plot(0:1/n150:1-1/n150, psd150)
        title(sprintf('Averaged periodogram for 150s of Trial %.0f', i))
        xlabel('normalised frequency (2\pi rad/sample)')
        ylabel('Magnitude (Au)')
        xlim([0 0.2])

    end
end


