function P = pgm(x)
    N = length(x);
    
    freq = 0;
    n = 1:N;
    exponential = exp(-2*pi*1i*n/N);
    
%     for pp = 1:N
%         P(pp) = 1/N * abs(sum(x .* (exponential.^freq)))^2;
%         freq = freq + 1/N;
%     end
%     
    X = fft(x);
    P = 1/N .* abs(X).^2;
    
    assert(length(P) == N, 'The Periodogram and the input signal are not of the same length')
