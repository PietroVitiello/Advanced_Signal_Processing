%% Part 1
clear all;
N = 10000;
v = randn(1, N);
[frequencies, centers, probs] = pdf(v, 30);

figure(), hold on
xlabel('values')
ylabel('probability')
bar(centers, probs);
plot(-5:0.1:5, normpdf(-5:0.1:5, 0, 1), 'r', 'LineWidth', 2);
legend('pdf()', 'theoretical pdf')

%% Part 2
addpath(genpath([fileparts(pwd), '\Subpart2']));

n_bins = 10;
syms x;
theoretical = piecewise(-1>x, 0, -1<x<2, 1/3, x>2, 0);

for i = 1:3
    [~, v3{i+3}, v3{i}] = pdf(rp3(1, 10^(i+1)), n_bins);
end

figure(2)
for i = 1:3
    subplot_title = sprintf('estimated pdf for %.0f samples', 10^(i+1));
    subplot(3, 1, i), hold on;
    bar(v3{i+3}, v3{i});
    fplot(theoretical,[-1.5 2.5], 'LineWidth', 2);
    title(subplot_title)
    xlabel('values');
    ylabel('probability')
end

%% Part 3
clear all;

non_st = [normrnd(0,1, 1, 50000) normrnd(10, 1, 1, 50000)];
time_axis = 1:100000;
[~, centers, probs] = pdf(non_st, 100);

figure()
plot(time_axis, non_st);
xlabel('Time samples (n)')
ylabel('Amplitude(Au)')

figure(), hold on
xlabel('values')
ylabel('probability')
bar(centers, probs);
