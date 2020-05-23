function [y, error, coefficients] = lms_signed_regressor(x, step, order)
    a = zeros(order, 1);
    for n = order+1:length(x)
        
        predictor = a' * x(n-1 : -1 : n-order);
        e = x(n) - predictor;
        a = a + step*e*sign(x(n-1 : -1 : n-order));
        
        y(n-order) = predictor;
        error(n-order) = e;
        coefficients(:,n-order) = a;
    end
end