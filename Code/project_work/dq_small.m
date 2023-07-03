clear
close all

% Parameters
R = 15;          %Resistance
L = 27e-3;       %Inductance
f_c = 50;       %Grid freq (Hz)
f_p = 10;       %Measurement freq (Hz)
tStep_min = 5e-6;   %Sim min step size
n = 10;         %No of measurement points

%%
f = sym('f');
Z_ana = [R+1i*2*pi*f*L, -2*pi*f_c*L
         2*pi*f_c*L, R+1i*2*pi*f*L];
s = tf('s');
Z_sys = [R+s*L, -2*pi*f_c*L
         2*pi*f_c*L, R+s*L];

freq_vals = logspace(2,5,50);
Z_dd = eval(subs(Z_ana(1,1),f,freq_vals));
Z_dq = eval(subs(Z_ana(1,2),f,freq_vals));
Z_qd = eval(subs(Z_ana(2,1),f,freq_vals));
Z_qq = eval(subs(Z_ana(2,2),f,freq_vals));

figure
subplot(4,2,1)
semilogx(freq_vals,20.*log10(abs(Z_dd)))
title('Z_{dd}')
grid on
subplot(4,2,3)
semilogx(freq_vals,angle(Z_dd)./(2*pi).*360)
grid on
subplot(4,2,2)
semilogx(freq_vals,20.*log10(abs(Z_dq)))
title('Z_{dq}')
grid on
subplot(4,2,4)
semilogx(freq_vals,angle(Z_dq)./(2*pi).*360)
grid on
subplot(4,2,5)
semilogx(freq_vals,20.*log10(abs(Z_qd)))
title('Z_{qd}')
grid on
subplot(4,2,7)
semilogx(freq_vals,rad2deg(angle(Z_qd)))
grid on
subplot(4,2,6)
semilogx(freq_vals,20.*log10(abs(Z_qq)))
title('Z_{qq}')
grid on
subplot(4,2,8)
semilogx(freq_vals,rad2deg(angle(Z_qq)))
grid on

%% Run model and gather data

tStep_temp = 1/f_p;
tStep = min(tStep_temp,tStep_min);

tStop = 1/f_p*10;
freq = [f_p f_p];
amp = [1, 0];
out1 = sim('dq0_fun3.slx');
amp = [0, 1];
out2 = sim('dq0_fun3.slx');

%% Plot data

Vd1 = 1;
Vq1 = 0;
Vd2 = 0;
Vq2 = 1;
Id1 = out1.I_d_mag(end,1)*exp(2i*pi*out1.I_d_phase(end,1)/360);
Id2 = out2.I_d_mag(end,1)*exp(2i*pi*out2.I_d_phase(end,1)/360);
Iq1 = out1.I_q_mag(end,1)*exp(2i*pi*out1.I_q_phase(end,1)/360);
Iq2 = out2.I_q_mag(end,1)*exp(2i*pi*out2.I_q_phase(end,1)/360);

vs = [Vd1 Vd2
      Vq1 Vq2];
is = [Id1 Id2
      Iq1 Iq2];
Z_meas = inv(is)*vs
Z = eval(subs(Z_ana,f,f_p))
Ydd = Id1/Vd1;
Yqd = Iq1/Vd1;
Ydq = Id2/Vq2;
Yqq = Iq2/Vq2;
Y_meas = [Ydd Yqd; Ydq Yqq]
Y = inv(eval(subs(Z_ana,f,f_p)))