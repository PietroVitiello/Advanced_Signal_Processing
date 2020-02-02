function [frequencies, centers, probs] = pdf(x, n_bins)
    N = length(x); %length of the realization
    if (N >= 100)
        [counts, centers] = hist(x, n_bins);
        frequencies = counts/N; %the percentage of points of the realization that fall in
                                %the range represented by the bin
        bin_width = centers(2) - centers(1); %how many values are represented by each bin
        probs = frequencies./bin_width; %represents the pdf, where each bin has been normalized
    else
        disp("Number of samples must be greater than 100")
    end
    
    area = sum(probs.*bin_width);
    if area == 1
        disp('True')
    else
        disp('The pdf does not integrate to 1')
        disp(bin_width)
    end
    

