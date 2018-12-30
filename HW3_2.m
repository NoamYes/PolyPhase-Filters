%%

clear all;
clc;
close all;

% this is a code for wet 3 of DSP

Fs = 8000; % original FS
Fs_y = 72000;

%% Section 1
t = linspace(0,pi,1000);
f1 = 1000;
f2 = 2000;
f3 = 3000;
phi1 = rand(1)*2*pi;
phi2 = rand(1)*2*pi;
phi3 = rand(1)*2*pi;
x = cos(2*pi*f1*t + phi1) + cos(2*pi*f2*t + phi2) + cos(2*pi*f3*t + phi3);

Ts = 1/Fs;
t_new = linspace(0, 511*Ts, 512);
x_n = cos(2*pi*f1*t_new + phi1) + cos(2*pi*f2*t_new + phi2) + cos(2*pi*f3*t_new + phi3);

Ts_y = 1/Fs_y;
t_new_y = linspace(0, 511*Ts_y, 512);
y_m = cos(2*pi*f1*t_new_y + phi1) + cos(2*pi*f2*t_new_y + phi2) + cos(2*pi*f3*t_new_y + phi3);

figure(1);
plot(t_new, x_n,'-', 'LineWidth', 1.5);
