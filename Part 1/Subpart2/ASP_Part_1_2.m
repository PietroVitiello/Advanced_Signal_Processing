% 
% n = 1:N;
% x = sin(2*pi*5*0.001*n)';
% eta = randn(length(n),1);
% y = x + eta;

%% Set up
clear all;
N = 100; %number of discrete times
M = 100; %number of realizations
x = rp1(M, N);
y = rp2(M, N);
v = rp3(M, N);

%Plot a realization
figure(), plot(1:N, x(1,:));

for i = 1: N
    ensamble_mx(i) = mean(x(:, i));
    ensamble_sx(i) = std(x(:, i));
    
    ensamble_my(i) = mean(y(:, i));
    ensamble_sy(i) = std(y(:, i));
    
    ensamble_mv(i) = mean(v(:, i));
    ensamble_sv(i) = std(v(:, i));
end

figure()
subplot(2, 1, 1)
plot(1:N, ensamble_mx)
subplot(2, 1, 2)
plot(1:N, ensamble_sx)

figure()
subplot(2, 1, 1)
plot(1:N, ensamble_my)
subplot(2, 1, 2)
plot(1:N, ensamble_sy)

figure()
subplot(2, 1, 1)
plot(1:N, ensamble_mv)
subplot(2, 1, 2)
plot(1:N, ensamble_sv)


