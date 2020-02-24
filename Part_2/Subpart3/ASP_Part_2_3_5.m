close all;
clear all;

load sunspot.dat
sun = sunspot(:,2);
sun = zscore(sun);