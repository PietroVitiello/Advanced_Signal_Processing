%% Set up
clear all;
N = 100; %number of discrete times
M = 100; %number of realizations
x = rp1(M, N);
y = rp2(M, N);
v = rp3(M, N);

%Plot a realization
%figure(), plot(1:N, x(1,:));

%%% Part 1
for i = 1: N
    ensemble_mx(i) = mean(x(:, i));
    ensemble_sx(i) = std(x(:, i));
    
    ensemble_my(i) = mean(y(:, i));
    ensemble_sy(i) = std(y(:, i));
    
    ensemble_mv(i) = mean(v(:, i));
    ensemble_sv(i) = std(v(:, i));
end

mean_estimatorx = mean(ensemble_mx);
std_estimatorx = mean(ensemble_sx);
mean_estimatory = mean(ensemble_my);
std_estimatory = mean(ensemble_sy);
mean_estimatorv = mean(ensemble_mv);
std_estimatorv = mean(ensemble_sv);

%% Plot Part 1
%Not stationary
figure()
subplot(2, 1, 1)
plot(1:N, ensemble_mx)
title('rp1')
xlabel('time samples')
ylabel('Ensemble mean')
subplot(2, 1, 2)
plot(1:N, ensemble_sx)
xlabel('time samples')
ylabel('Ensemble standard deviation')

%stationary
figure()
subplot(2, 1, 1)
plot(1:N, ensemble_my)
title('rp2')
xlabel('time samples')
ylabel('Ensemble mean')
subplot(2, 1, 2)
plot(1:N, ensemble_sy)
xlabel('time samples')
ylabel('Ensemble standard deviation')

%stationary
figure()
subplot(2, 1, 1)
plot(1:N, ensemble_mv)
title('rp3')
xlabel('time samples')
ylabel('Ensemble mean')
subplot(2, 1, 2)
plot(1:N, ensemble_sv)
xlabel('time samples')
ylabel('Ensemble standard deviation')

%% Part 2
%clear all;
N = 1000; %number of discrete times
M = 4; %number of realizations
x = rp1(M, N);
y = rp2(M, N);
v = rp3(M, N);

for i = 1: M
    %not ergodic
    sample_mx(i, 1) = mean(x(i, :));
    sample_sx(i, 1) = std(x(i, :));
    
    %not ergodic
    sample_my(i, 1) = mean(y(i, :));
    sample_sy(i, 1) = std(y(i, :));
    
    %ergodic
    sample_mv(i, 1) = mean(v(i, :));
    sample_sv(i, 1) = std(v(i, :));
end





