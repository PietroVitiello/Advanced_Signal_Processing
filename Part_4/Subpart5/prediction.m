function predicted = prediction(x, a, do_figure)
    order = length(a);
    
    for n = 1:length(x)-order    
        predicted(n) = a' * x(n+order-1 : -1 : n);   
    end
    
    if do_figure == 1
        figure(), hold on
        plot(predicted)
        plot(x(order+1:end))
    end
end