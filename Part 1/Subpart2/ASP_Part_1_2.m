%% Set up
clear all;
N = 100; %number of discrete times
M = 100; %number of realizations
x = rp1(M, N);
y = rp2(M, N);
v = rp3(M, N);

%Plot a realization
%figure(), plot(1:N, x(1,:));

for i = 1: N
    ensamble_mx(i) = mean(x(:, i));
    ensamble_sx(i) = std(x(:, i));
    
    ensamble_my(i) = mean(y(:, i));
    ensamble_sy(i) = std(y(:, i));
    
    ensamble_mv(i) = mean(v(:, i));
    ensamble_sv(i) = std(v(:, i));
end

mean_estimatorx = mean(ensamble_mx);
std_estimatorx = mean(ensamble_sx);
mean_estimatory = mean(ensamble_my);
std_estimatory = mean(ensamble_sy);
mean_estimatorv = mean(ensamble_mv);
std_estimatorv = mean(ensamble_sv);

%Not stationary
figure()
subplot(2, 1, 1)
plot(1:N, ensamble_mx)
subplot(2, 1, 2)
plot(1:N, ensamble_sx)

%stationary
figure()
subplot(2, 1, 1)
plot(1:N, ensamble_my)
subplot(2, 1, 2)
plot(1:N, ensamble_sy)

%stationary
figure()
subplot(2, 1, 1)
plot(1:N, ensamble_mv)
subplot(2, 1, 2)
plot(1:N, ensamble_sv)

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

%% Part 3
N = 100; %number of discrete times
M = 100; %number of realizations
x_math = rp1_mathematical(M, N);

figure()
plot(1:N, ensamble_mx, x_math)




