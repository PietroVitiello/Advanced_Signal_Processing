function [out, error, weight_matrix] = lms(x, z, adaptation_gain, order)
    N_w = order - 1;
    w = zeros(order, 1);
    step = adaptation_gain;
    
    iterations = length(x) - N_w;
    for i = 1:iterations
        weight_matrix(:,i) = w;
        
        y = w' * x(i+N_w : -1 : i);
        e = z(i+N_w) - y;
        w = w + step*e*x(i+N_w : -1 : i);
        
        error(i) = e;
        out(i) = y;
        
    end
    
    weight_matrix(:,i) = w;
end