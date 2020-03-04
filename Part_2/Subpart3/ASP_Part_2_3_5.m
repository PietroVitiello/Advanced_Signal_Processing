%% Prediction
close all;
clear all;

N = 140;
axis = 1:N;

load sunspot.dat
sun = sunspot(axis,2);
sun = zscore(sun);

prediction = {};
i = 1;
for p = [1 2 10]
    model = ar(sun, p);
    ii = 1;
    for h = [1 2 5 10]
        prediction{i}(ii, :) = predict(model, sun, h);
        
        figure(), hold on;
        plot(axis, prediction{i}(ii, :), 'r')
        plot(axis, sun, 'k')
        title(sprintf("Predicted sunspot data for order %.0f and horizon %.0f", p, h));
        xlabel("sample (n)")
        ylabel("amplitude (Au)")
        
        ii = ii + 1;
    end
    i = i + 1;
end

%% Model comparison
a = aryule(sun, 2);
noise = randn(140, 1);
filtered_sun = filter(-a, 1, sun);
axis = 0:length(sun)-1;

model = ar(sun, 2);
predictionn = predict(model, sun, 0);

figure(), hold on;
plot(axis, sun)
plot(axis, predictionn)




