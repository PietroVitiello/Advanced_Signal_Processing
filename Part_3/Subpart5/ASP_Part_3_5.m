clear all
close all
clc

%% Set up
pwdd = erase(fileparts(pwd), 'Part_3');
addpath(genpath([pwdd, 'Part_2\Subpart5']));

load RRI_2_5.mat

which = [0];

% plot(RRI_trial3)

trials = {RRI_trial1; RRI_trial2; RRI_trial3};
frequencies = [RRI_fs1 RRI_fs2 RRI_fs3];
clear RRI_trial1 RRI_trial2 RRI_trial3 RRI_fs1 RRI_fs2 RRI_fs3

if (which >= 1 && which <=3)
    for i = which
        [psd, psd50, psd150, n50, n150] = psd_trials(trials{i}, frequencies(i));

        figure()
        plot(0:1/length(psd):1-1/length(psd), psd)
        title(sprintf('Standard periodogram for Trial %.0f', i))
        xlim([0 0.1])

        figure()
        plot(0:1/n50:1-1/n50, psd50)
        title(sprintf('Averaged periodogram for 50s of Trial %.0f', i))
        xlim([0 0.1])

        figure()
        plot(0:1/n150:1-1/n150, psd150)
        title(sprintf('Averaged periodogram for 150s of Trial %.0f', i))
        xlim([0 0.1])

    end
end

% figure(1)
% plot((0:length(standard_psd)/2 - 1)/length(standard_psd),standard_psd(1:end/2))
% ylim([0 40]) % For all trials
% % title(['PSD of RRI Trial ' num2str(trial_num)])
% xlabel('Frequency')
% ylabel('PSD')
% 
% figure(2)
% plot((0:99)/200, psd50(1:100))
% % title(['PSD of RRI Trial ' num2str(trial_num) ' (50 seconds average)'])
% xlabel('Normalized Frequency')
% ylabel('PSD')
% grid on
% grid minor
% 
% figure(3)
% plot((0:299)/600, psd150(1:300))
% % title(['PSD of RRI Trial ' num2str(trial_num) ' (150 seconds average)']);
% xlabel('Normalized Frequency')
% ylabel('PSD')
% grid on
% grid minor


