%% Part 1
clear all;
N = 1000;
v = randn(1, N);
[frequencies, centers] = pdf(v, 20);
bar(centers, frequencies);

%% Part 2
addpath(genpath([fileparts(pwd), '\Subpart2']));

n_bins = 20;
syms x;
theoretical = piecewise(-1>x, 0, -1<x<2, 1/n_bins, x>2, 0);

for i = 1:3
    [v3{i}, v3{i+3}] = pdf(rp3(1, 10^(i+1)), n_bins);
end

figure()
for i = 1:4
    if i == 4
        subplot(4, 1, 4)
        fplot(theoretical,[-1.5 2.5])
        ylim([0 1.2/n_bins])
    else
        subplot(4, 1, i)
        bar(v3{i+3}, v3{i}); 
    end
end

%% Part 3

    