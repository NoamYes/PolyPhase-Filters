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


theta = 0:0.001:1;
transition = fstop - fpass;
line = (-1/transition)+fstop/transition;
lpFilt = (theta < fpass)*L + L*(theta > fpass & theta < fstop).*(-1/transition*theta+fstop/transition);
figure(1)
plot(theta, (lpFilt),'-', 'LineWidth', 1.5);
xlim([ 0 1 ])
title('Ideal Prototype lowpass interpolation filter') 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response')


%% Section 2 

% equiripple remez
N = 100;
b = firpm(N-1, theta(1:end-1), lpFilt(1:end-1));
[h,w] = freqz(b,1,512);
figure(2)
plot(theta,lpFilt, w/pi,abs(h));
xlim([ 0 1 ])
title(['Equiripple Remez lowpass polyphase filter, Order = ' num2str(N)]) 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response')

Phases_Mat_coeffs = zeros(L,N);

figure(3)

for i=0:L-1
    %Phases_Mat_coeffs(i+1,:) = circshift(b, i);
    shifted = circshift(b, i);
    grpdelay(shifted(1:L:end))
    hold on;
end

ylim([-10 20]);
title(['Phase Delays of polyphase filters, Order = ' num2str(N)])

figure(4) 
legend = {};
for i=0:L-1
    %Phases_Mat_coeffs(i+1,:) = circshift(b, i);
    shifted = circshift(b, i);
    [h,w] = freqz(shifted(1:L:end),1);
    plot(w/pi,20*log10(abs(h)))
    %legappend(['Pol ',num2str(i)])
    hold on;
end


xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
