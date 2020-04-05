function [out, error, weight_matrix] = lms_gear(x, z, adaptation_gain, scaling, threshold , order)
    N_w = order - 1;
    w = zeros(order, 1);
    step = adaptation_gain;
    avg_mse = 0;
    
    iterations = length(x) - N_w;
    for i = 1:iterations
        weight_matrix(:,i) = w;
        
        y = w' * x(i+N_w : -1 : i);
        e = z(i+N_w) - y;
        
        if i<=5
            avg_mse = avg_mse + (e^2)/5;
            reference_mse = avg_mse;
        else
            avg_mse = avg_mse - (error(i-5))^2/5 + (e^2)/5;
        end
        
        if (avg_mse < reference_mse/threshold)
            step = step * (1-scaling/100);
            reference_mse = avg_mse;
        end
        
        w = w + step*e*x(i+N_w : -1 : i);
        
        error(i) = e;
        out(i) = y;
        
    end
    
    weight_matrix(:,i) = w;
end