function [y, error, coefficients] = speech_lms(x, z, step, order)
    a = zeros(order, 1);
    for n = order+1:length(x)
        
        predictor = a' * x(n-1 : -1 : n-order);
        e = z(n) - predictor;
        a = a + step*e*x(n-1 : -1 : n-order);
        
        y(n-order) = predictor;
        error(n-order) = e;
        coefficients(:,n-order) = a;
    end
end