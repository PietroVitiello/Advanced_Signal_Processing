function [frequencies, centers] = pdf(x, n_bins)
    N = length(x);
    if (N >= 100)
        [counts, centers] = hist(x, n_bins);
        frequencies = counts/N;
    else
        disp("Number of samples must be greater than 100")
    end
    

