%% set up
clear all;
a1 = rand(100, 1).*5 - 2.5;
a2 = rand(100, 1).*3 - 1.5;
w = randn(1000, 1);

%% Part 1
clear pairs
tol = 10000;
pair_idx = 1;

for i = 1:length(a1)
    for j = 1:length(a2) 
        a = [1 -a1(i) -a2(j)];
        x = filter(1, a, w);
        
        if (x(end) < tol) && (x(end) > -tol)
            pairs(pair_idx,:) = [a1(i), a2(j)];
            pair_idx = pair_idx+1;
        end
    end 
end

scatter(pairs(:,1), pairs(:,2), 'r*')
title("Stable coefficients pairs for 100 coefficient values")
xlabel("a_{1}")
ylabel("a_{2}")
xlim([-2.5 2.5])
ylim([-1.5 1.5])


