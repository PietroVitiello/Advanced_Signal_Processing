%% Creation of signal
clc, clear all, close all;
n_samples = 1000;
which_exercise = 1; %choose 1 for uniform distribution and 2 for gaussian distribution
x_axis = 1:n_samples;

if (which_exercise == 1)
    X = rand(n_samples, 1);  %Uniformly distributed between 0 and 1
elseif (which_exercise == 2)
    X = randn(n_samples, 1);  %Gaussian distribution with 0 mean and unit standard deviation
end

%% Plot the signal
figure()
plot(x_axis, X)
ylabel("Amplitude (Au)");
xlabel("Samples (n)");

%% Part 1
sample_mean = mean(X);

%% Part 2
sample_standard_deviation = std(X);

%% Part 3
realizations = 1:10;

if (which_exercise == 1)
    x_realizations = rand(10, n_samples);
elseif (which_exercise == 2)
    x_realizations = randn(10, n_samples);
end

sample_mean_realizations = zeros(10, 1);
std_realizations = zeros(10, 1);
error_realizations = zeros(10, 1);

for realization = realizations  
    sample_mean_realizations(realization, 1) = mean(x_realizations(realization, :));
    std_realizations(realization, 1) = std(x_realizations(realization, :));
    if (which_exercise == 1)
        error_realizations(realization, 1) = 0.5 - sample_mean_realizations(realization, 1);
    elseif (which_exercise == 2)
        error_realizations(realization, 1) = 0 - sample_mean_realizations(realization, 1);
    end
end

m_estimatorBias = 0.5 - mean(sample_mean_realizations);
s_estimatorBias = 1/sqrt(12) - mean(std_realizations);

%% Part 3 plotting
figure()

subplot(2,1,1)
hold on
scatter(realizations, sample_mean_realizations, 18,'bo', 'filled')
title("realizations' sample mean")
ylabel("sample mean")
xlabel("realizations")
%ylim([0.45 0.55])
if (which_exercise == 1)
    hline = refline([0, 0.5]);
    hline.Color = 'b';
    for i = realizations
        plot([i, i], [0.5, sample_mean_realizations(i)], 'b-')
    end
elseif (which_exercise == 2)
    hline = refline([0, 0]);
    hline.Color = 'b';
    for i = realizations
        plot([i, i], [0, sample_mean_realizations(i)], 'b-')
    end
end

subplot(2,1,2)
hold on
scatter(realizations, std_realizations, 18, 'bo', 'filled')
title("realizations' sample standard deviation")
ylabel("sample standard deviation")
xlabel("realizations")
%ylim([0.25 0.35])
if (which_exercise == 1) 
    hline = refline([0, 1/sqrt(12)]);
    hline.Color = 'b';
    for i = realizations
        plot([i, i], [1/sqrt(12), std_realizations(i)], 'b-')
    end
elseif (which_exercise == 2)
    hline = refline([0, 1]);
    hline.Color = 'b';
    for i = realizations
        plot([i, i], [1, std_realizations(i)], 'b-')
    end
end

%% Part 4
n_bins = 3; %n_samples;
title_text = sprintf("histogram pdf for %.0f samples and %.0f bins", n_samples, n_bins);
title_histogram = sprintf("pdf estimate using the 'Normalized' histogram() function");
figure()
hold on
title(title_text);
ylabel('Probability')
[counts, centers] = hist(X, n_bins);
bar(centers, counts/n_samples)

if (which_exercise == 1)
    hline = refline([0, 1/n_bins]);
    hline.Color ='r';
    xlabel('x');
elseif (which_exercise == 2)
    xlabel('y');
    figure(), hold on
    histogram(X, 'Normalization', 'pdf');
    plot(-5:0.1:5, normpdf(-5:0.1:5, 0, 1));
    title(title_histogram);
    xlabel('y');
    ylabel('Probability')
end





