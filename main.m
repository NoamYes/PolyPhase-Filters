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


%theta = 0:0.001:1;
transition = fstop - fpass;
line = (-1/transition)+fstop/transition;
%lpFilt = (theta < fpass)*L + L*(theta > fpass & theta < fstop).*(-1/transition*theta+fstop/transition);
lpFilt = [9 9 0 0];
theta = [0 0.7/L 1.3/L 1];

figure(1)
plot(theta, (lpFilt),'-', 'LineWidth', 1.5);
xlim([ 0 1 ])
title('Ideal Prototype lowpass interpolation filter') 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response')


%% Section 2 

% equiripple remez
N = 135;
b = firpm(N-1, theta, lpFilt);
[h,w] = freqz(b,1,512);
figure(2)
% plot(theta,lpFilt, w/pi,10*log10(abs(h)));
plot(w/pi,20*log10(abs(h)));
h_ripple = h;
hold on;
xlim([ 0 1 ])
title(['Equiripple Remez lowpass polyphase filter, Order = ' num2str(N)]) 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response[dB]')
legend({'Equirriple'})

polyPhaseEquiripple=cell(1,L);

figure(3)

legendPol = {};
for i=0:L-1
    %Phases_Mat_coeffs(i+1,:) = circshift(b, i);
    shifted = circshift(b, i);
    polyPhaseEquiripple{i+1} = shifted(1:L:end);
    grpdelay(shifted(1:L:end))
    legendPol{end+1} = ['Pol ' num2str(i)];
    hold on;
end
legend(legendPol);
ylim([-10 20]);
title(['Phase Delays of polyphase filters, Equiripple, Order = ' num2str(N)])

figure(4) 
legendPol = {};
for i=0:L-1
    %Phases_Mat_coeffs(i+1,:) = circshift(b, i);
    shifted = circshift(b, i);
    [h,w] = freqz(shifted(1:L:end),1);
    plot(w/pi,20*log10(abs(h)))
    legendPol{end+1} = ['Pol ' num2str(i)];
    hold on;
end
legend(legendPol);
title(['Magnitude Response of polyphase filters, Equiripple, Order = ' num2str(N)])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% Section 3

% LS

N = 100;
b = firls(N-1, theta, lpFilt);
[h,w] = freqz(b,1,512);
h_ls = h;
figure(5)
plot(w/pi,10*log10(abs(h_ls)));
hold on;
plot(w/pi,10*log10(abs(h_ripple)));
xlim([ 0 1 ])
title(['LS and Equirriple polyphase filter, Order = ' num2str(N)]) 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response')
legend({'LS', 'Equiripple'});

polyPhaseLS = cell(1,L);

figure(6)

legendPol = {};
for i=0:L-1
    %Phases_Mat_coeffs(i+1,:) = circshift(b, i);
    shifted = circshift(b, i);
    polyPhaseLS{i+1} = shifted(1:L:end);
    grpdelay(shifted(1:L:end))
    legendPol{end+1} = ['Pol ' num2str(i)];
    hold on;
end
legend(legendPol);
ylim([-10 20]);
title(['Phase Delays of polyphase filters, LS, Order = ' num2str(N)])
figure(7) 
legendPol = {};
for i=0:L-1
    %Phases_Mat_coeffs(i+1,:) = circshift(b, i);
    shifted = circshift(b, i);
    [h,w] = freqz(shifted(1:L:end),1);
    plot(w/pi,20*log10(abs(h)))
    legendPol{end+1} = ['Pol ' num2str(i)];
    hold on;
end
legend(legendPol);
title(['Magnitude Response of polyphase filters, LS, Order = ' num2str(N)])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')


%% this is a code for wet 3 Q2 of DSP

Fs = 8000; % original FS
Fs_y = 72000;

%% Q2 Section 1
t = linspace(0,pi,1000);
f1 = 1000;
f2 = 2000;
f3 = 3000;
phi1 = unifrnd(0, 2*pi);
phi2 = unifrnd(0, 2*pi);
phi3 = unifrnd(0, 2*pi);

Ts = 1/Fs;
t_new = 0:Ts:(511*Ts);
x_n = cos(2*pi*f1*t_new + phi1) + cos(2*pi*f2*t_new + phi2) + cos(2*pi*f3*t_new + phi3);

Ts_y = 1/Fs_y;
t_new_y = 0:Ts_y:(511*Ts_y);
y_m = cos(2*pi*f1*t_new_y + phi1) + cos(2*pi*f2*t_new_y + phi2) + cos(2*pi*f3*t_new_y + phi3);

figure(8);
plot(t_new, x_n,'-', 'LineWidth', 1.5);
%plot(t, x,'-', 'LineWidth', 1.5);
hold on;
%scatter(t_new, x_n);
scatter(t_new_y, y_m);
title('x(t) sampled in 8khz and 72khz');
xlim([0, 5e-3])
legend('x_n 8khz', 'y_m 72khz');

%% Q2  Section 2

% EQUIRIPPLE RECOVERY 

u = cell(1,L);
y_length = length(upsample(filter(polyPhaseEquiripple{1},1,x_n),L));
y_hat = zeros(1,y_length);
for branch = 0:L-1 
    
    x_conv_p = filter(polyPhaseEquiripple{branch+1},1,x_n);
    ui = upsample(x_conv_p, L);
    y_hat_i = circshift(ui, branch);
    y_hat_i(1:branch) = 0;
    y_hat = y_hat + y_hat_i;
    
end

figure(9);

plot(y_m, '-', 'LineWidth', 1.5);
hold on;
plot(y_hat, '-', 'LineWidth', 1.5);
xlim([0,512]);


%% Q2 Section 3
% EQUIRIPPLE RECOVERY 

u = cell(1,L);
y_length = length(upsample(filter(polyPhaseLS{1},1,x_n),L));
y_hat_ls = zeros(1,y_length);
for branch = 0:L-1 
    
    x_conv_p = filter(polyPhaseLS{branch+1},1,x_n);
    ui = upsample(x_conv_p, L);
    y_hat_i = circshift(ui, branch);
    y_hat_i(1:branch) = 0;
    y_hat_ls = y_hat_ls + y_hat_i;
    
end

figure(10)

plot(y_m, '-', 'LineWidth', 1.5);
hold on;
plot(y_hat_ls, '-', 'LineWidth', 1.5);
xlim([0,512]);
