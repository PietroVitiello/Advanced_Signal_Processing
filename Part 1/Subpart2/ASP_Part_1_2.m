% N = 1000; %number of discrete times
% M = 10 %
% 
% n = 1:N;
% x = sin(2*pi*5*0.001*n)';
% eta = randn(length(n),1);
% y = x + eta;

%% Set up
clear all;
x = rp1(100, 100);
y = rp2(100, 100);
v = rp3(100, 100);