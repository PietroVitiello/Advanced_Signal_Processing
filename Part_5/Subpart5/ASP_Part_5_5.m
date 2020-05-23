clear all
close all
clc

pwdd = erase(fileparts(pwd), 'Part_5');
addpath(genpath([pwdd, 'Part_3\Subpart1']));

%% Creating the signal
n_frequencies = 500;
frequencies = linspace(0, 0.5, n_frequencies);
estimated_f = zeros(1, n_frequencies);

N = 10;
n = 0 : N-1;
do_per_figures = 0;

figure(), hold on
for f = frequencies
    x = cos(2*pi*f*n);
    
    per = pgm(x);
    norm_freq = linspace(0, 0.5, length(per)/2);
    [~, estimated_f(find(frequencies == f))] = max(per(1:end/2));
    estimated_f(find(frequencies == f)) = norm_freq(estimated_f(find(frequencies == f))); 
    
    if (do_per_figures == 1)
        plot(norm_freq, per(1:end/2))
        title(['frequency ' num2str(f)])
    end
end

figure(), hold on
scatter(1:n_frequencies, frequencies, 'ob')
scatter(1:n_frequencies, estimated_f, '*r')


