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


% theta = 0:0.001:1;
% transition = fstop - fpass;
% line = (-1/transition)+fstop/transition;
% lpFilt = (theta < fpass)*L + L*(theta > fpass & theta < fstop).*(-1/transition*theta+fstop/transition);

theta = [0 fpass fstop 1];
lpFilt = [9 9 0 0];

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
figure(2);
plot(w/pi,db(abs(h)));
h_ripple = h;
hold on;
xlim([ 0 1 ])
title(['Equiripple Remez lowpass polyphase filter, Order = ' num2str(N)]) 
xlabel('\theta [\pi rads]')
ylabel('Amplitude Response')
legend({'Equirriple'})

polyPhaseEquiripple=cell(1,L);

figure(3);

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

N = 135;
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






% this is a code for wet 3 Q2 of DSP

Fs = 8000; % original FS
Fs_y = 72000;

%% Q2 Section 1
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
t_new_y = 0:Ts_y:511*Ts;
y_m = cos(2*pi*f1*t_new_y + phi1) + cos(2*pi*f2*t_new_y + phi2) + cos(2*pi*f3*t_new_y + phi3);

figure(8);
plot(t_new, x_n,'-', 'LineWidth', 1.5);
hold on;
scatter(t_new, x_n);
scatter(t_new_y, y_m);
xlim([0, 5e-3]);
title('x[n] and y[m]');
legend({'x[n]', 'y[m]'});


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
xlim([1 512])
title('y[m] and $\hat{y}[m]$','Interpreter','latex');
legend({'y[m]','$\hat{y}[m]$'},'Interpreter','latex');

corrEq = xcorr(y_hat, y_m);
[~,tmp] = max(corrEq);
delay_Eq = mod(tmp,y_length);
estD = delay_Eq;

[minY_m,~] = min(y_m);
[minY_eq,~] = min(y_hat);
[maxY_m,~] = max(y_m);
[maxY_eq,~] = max(y_hat);

shifted_eq = circshift(y_m, estD);
tmp = shifted_eq(1:length(y_m));
diffY_Eq = y_m - tmp;
figure(10)
plot(diffY_Eq, '-', 'LineWidth', 1.5);
title('Equiripple difference y[m] and $\hat{y}[m]$','Interpreter','latex');
SE_eq = sum((tmp - y_m).^2);


%% Q2 Section 3
% LS RECOVERY 

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

figure(11)

plot(y_m, '-', 'LineWidth', 1.5);
hold on;
plot(y_hat_ls, '-', 'LineWidth', 1.5);
title('y[m] and $\hat{y}[m]$','Interpreter','latex');
legend({'y[m]','$\hat{y}[m]$'},'Interpreter','latex');
xlim([1 512])

corrLS = xcorr(y_hat_ls, y_m);
[~,tmp] = max(corrLS);
delay_Eq = mod(tmp,y_length);
estD = delay_Eq;
[minY_m,~] = min(y_m);
[minY_LS,~] = min(y_hat_ls);
[maxY_m,~] = max(y_m);
[maxY_LS,~] = max(y_hat_ls);

shifted_eq = circshift(y_m, estD);
tmp = shifted_eq(1:length(y_m));
diffY_LS = y_m - tmp;
figure(12)
plot(diffY_LS, '-', 'LineWidth', 1.5);
title('LS difference y[m] and $\hat{y}[m]$','Interpreter','latex');
SE_LS = sum((tmp - y_m).^2);


