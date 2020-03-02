clc
clear all
close all

load RAW.mat

%show data
figure()
plot(data)

%find breaking points
[~, first_clap] = min(data(240000:260000));
first_clap = 240000 + first_clap;
[~, second_clap] = min(data(500000:510000));
second_clap = 500000 + second_clap;

%divide the three trials
ECG_trial1 = data(1:first_clap-1);
ECG_trial2 = data(first_clap+1:second_clap-1);
ECG_trial3 = data(second_clap+1:end);

%show the three trials
figure()
subplot(3,1,1)
plot(ECG_trial1)
subplot(3,1,2)
plot(ECG_trial2)
subplot(3,1,3)
plot(ECG_trial3)

%convert the three trials into RRI
[RRI_trial1, RRI_fs1] = ECG_to_RRI(ECG_trial1, fs);
[RRI_trial2, RRI_fs2] = ECG_to_RRI(ECG_trial1, fs);
[RRI_trial3, RRI_fs3] = ECG_to_RRI(ECG_trial1, fs);

fprintf('\n\nConversion done')

%Save the trials
variables_to_save = {
    'RRI_trial1', 'RRI_fs1', 'RRI_trial2', 'RRI_fs2', 'RRI_trial3', 'RRI_fs3', 'data', 'fs'
};
save('RRI_2_5.mat', variables_to_save{:})

fprintf('Saved variables')
