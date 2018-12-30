%%

clear all;
clc;
close all;

% this is a code for wet 3 of DSP

Fs = 1000; % original FS

%% Section 1

% ideal filter prototype

L = 9;
fsamp = Fs*L;
M = 1000;
fstop = 1.3/L;
fpass = 0.7/L;


theta = 0:0.001*pi:pi;
transition = fstop - fpass;
line = (-1/transition)+fstop/transition;
lpFilt = (theta < fpass) + (theta > fpass & theta < fstop).*(-1/transition*theta+fstop/transition);
plot(theta, (lpFilt),'-', 'LineWidth', 1.5);
xlim([ 0 1 ])
title('Ideal Prototype lowpass interpolation filter') 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response')


%% Section 2 

% equiripple remez

