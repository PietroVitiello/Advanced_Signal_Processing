%% Creation of signal
clc; clear all;
n_samples = 1000;
which_exercise = 2; %choose 1 for uniform distribution and 2 for gaussian distribution
x_axis = 1:n_samples;

if (which_exercise == 1)
    X = rand(n_samples, 1);  %Uniformly distributed between 0 and 1
elseif (which_exercise == 2)
    X = randn(n_samples, 1);  %Gaussian distribution with 0 mean and unit standard deviation
end

figure()
plot(x_axis, X)
xlabel("discrete time (n)");

%% Part 1
sample_mean = mean(X);

%% Part 2
sample_standard_deviation = std(X);

%% Part 3
realizations = 1:10;
x_realizations = zeros(1000, 10);
sample_mean_realizations = zeros(1, 10);
std_realizations = zeros(1, 10);
bias_realizations = zeros(1, 10);

for realization = realizations
    if (which_exercise == 1)
        x_realizations(:, realization) = rand(1, 1000);
    elseif (which_exercise == 2)
        x_realizations(:, realization) = randn(n_samples, 1);
    end
    
    sample_mean_realizations(1, realization) = mean(x_realizations(:, realization));
    std_realizations(1, realization) = std(x_realizations(:, realization));
    bias_realizations(1, realization) = 0.5 - sample_mean_realizations(1, realization);
end


figure()
if (which_exercise == 1)
    subplot(2,1,1)
    scatter(realizations, sample_mean_realizations, 18,'ro', 'filled')
    title("realizations' sample means")
    ylabel("sample mean")
    xlabel("realization")
    ylim([0.45 0.55])
    hline = refline([0, 0.5]);

    subplot(2,1,2)
    scatter(realizations, std_realizations, 18, 'ro', 'filled')
    title("realizations' sample standard deviation")
    ylabel("sample standard deviation")
    xlabel("realization")
    ylim([0.25 0.35])
    hline = refline([0, 1/sqrt(12)]);
elseif (which_exercise == 2)
    subplot(2,1,1)
    scatter(realizations, sample_mean_realizations, 18,'ro', 'filled')
    title("realizations' sample means")
    ylabel("sample mean")
    xlabel("realization")
    ylim([-0.1 0.1])
    hline = refline([0, 0]);

    subplot(2,1,2)
    scatter(realizations, std_realizations, 18, 'ro', 'filled')
    title("realizations' sample standard deviation")
    ylabel("sample standard deviation")
    xlabel("realization")
    ylim([0.9 1.1])
    hline = refline([0, 1]);
end

%% Part 4
n_bins = 10;%n_samples;
figure()
hold on
[counts, centers] = hist(X, n_bins);
bar(centers, counts/n_samples)
normpdf((-n_samples:1/n_samples:n_samples), 0, 1);

if (which_exercise == 1)
    hline = refline([0, 1/n_samples]);
    hline.Color ='r';
elseif (which_exercise == 2)
    plot(-5:0.1:5, normpdf(-5:0.1:5, 0, 1));
end





