close all
T = 4e-5;
Fs = 1/T;            % Sampling frequency                    
L = 30000;             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

x1 = cos(2*pi*20*t-90*2*pi/360);
x2 = cos(2*pi*50*t+pi/4);
x3 = cos(2*pi*100*t);
X = x1 + x2 + x3;

a = 2*fft(x1)/L;
X_1 = interp1(f,a(1:L/2+1),20)
a = 2*fft(x2)/L;
X_2 = interp1(f,a(1:L/2+1),50)
a = 2*fft(x3)/L;
X_3 = interp1(f,a(1:L/2+1),100)

figure
plot(t,x1,t,x2,t,x3)
grid on
legend('x1','x2','x3')
xlim([0 0.05])

n_cycles = 50;
T_m = 1/100;
X2 = X(L-round((T_m*n_cycles)/T)+1:end);

figure
plot(t(L-round((T_m*n_cycles)/T)+1:end),X2)

L = length(X2);
Y = fft(X2);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 150])

Y2=Y;
threshold = max(abs(Y))/10000; %tolerance threshold
Y2(abs(Y)<threshold) = 0; %maskout values that are below the threshold

figure
plot(f,rad2deg(angle(Y2(1:L/2+1))))
xlim([0 150])
