clear
close all

% Parameters
R = 15;              %Resistance
L = 27e-3;           %Inductance
f_c = 50;           %Grid freq (Hz)
f_min = 18.853;         %Min measurement freq (Hz)
f_max = 10^6;       %Max measurement freq (Hz)
n = 20;             %No of measurement points
tStep_min = 5e-6;   %Sim min step size
tSettle = 50e-3;    %Settling time
amplitude = 100;    %Injection amplitude

%%
f = sym('f');
Z_ana = [R+1i*2*pi*f*L, -2*pi*f_c*L
         2*pi*f_c*L, R+1i*2*pi*f*L];
s = tf('s');
Z_sys = [R+s*L, -2*pi*f_c*L
         2*pi*f_c*L, R+s*L];

freq_vals = logspace(log10(f_min),log10(f_max),200);
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

freq_meas = logspace(log10(f_min),log10(f_max),n);
for i = 1:1:length(freq_meas)
    i
    %tStep_temp = 1/freq_meas(1,i);
    %tStep = min(tStep_temp,tStep_min);
    %tStep = 1/(freq_meas(1,i)*50);
    tStep = tStep_min;
    
    tStop = max(1/freq_meas(1,i)*10,tSettle);
    freq = [freq_meas(1,i) freq_meas(1,i)];
    amp = [amplitude, 0];
    out1 = sim('dq0_fun3.slx');
    amp = [0, amplitude];
    out2 = sim('dq0_fun3.slx');
    Vd1 = amplitude;
    Vq1 = 0;
    Vd2 = 0;
    Vq2 = amplitude;
    Id1 = out1.I_d_mag(end,1)*exp(2i*pi*out1.I_d_phase(end,1)/360);
    Id2 = out2.I_d_mag(end,1)*exp(2i*pi*out2.I_d_phase(end,1)/360);
    Iq1 = out1.I_q_mag(end,1)*exp(2i*pi*out1.I_q_phase(end,1)/360);
    Iq2 = out2.I_q_mag(end,1)*exp(2i*pi*out2.I_q_phase(end,1)/360);

    vs = [Vd1 Vd2
          Vq1 Vq2]
    is = [Id1 Id2
          Iq1 Iq2]
    Z_meas = inv(is)*vs
    Z_dd_meas(1,i) = Z_meas(1,1);
    Z_dq_meas(1,i) = Z_meas(1,2);
    Z_qd_meas(1,i) = Z_meas(2,1);
    Z_qq_meas(1,i) = Z_meas(2,2);
    Z = eval(subs(Z_ana,f,freq_meas(1,i)))
%     Z = eval(subs(Z_ana,f,f_p));
%     Ydd = Id1/Vd1;
%     Yqd = Iq1/Vd1;
%     Ydq = Id2/Vq2;
%     Yqq = Iq2/Vq2;
%     Y_meas = [Ydd Yqd; Ydq Yqq];
%     Y = inv(eval(subs(Z_ana,f,f_p)));
    clear out1 out2
end

%% Plot data
close all

figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_dd_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,20.*log10(abs(Z_dd)))
title('Z_{dd}')
ylabel('Magnitude (dB)')
legend('Measured','Analytical','Location','southeast')
xticklabels({})
ylim([20 120])
grid on
hold off
subplot(2,1,2)
semilogx(freq_meas,angle(Z_dd_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,angle(Z_dd)./(2*pi)*360)
legend('Measured','Analytical','Location','southeast')
ylabel('Phase (deg)')
xlabel('Frequency (Hz)')
grid on
hold off

Z_dq_meas(1,18) = conj(Z_dq_meas(1,18));
Z_dq_meas(1,19) = conj(Z_dq_meas(1,19));

figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_dq_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,20.*log10(abs(Z_dq)))
legend('Measured','Analytical','Location','southeast')
title('Z_{dq}')
ylabel('Magnitude (dB)')
xticklabels({})
ylim([18 19])
grid on
hold off
subplot(2,1,2)
semilogx(freq_meas,angle(Z_dq_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,angle(Z_dq)./(2*pi)*360)
legend('Measured','Analytical','Location','southeast')
ylabel('Phase (deg)')
xlabel('Frequency (Hz)')
ylim([175 185])
grid on
hold off

figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_qd_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,20.*log10(abs(Z_qd)))
legend('Measured','Analytical','Location','southeast')
title('Z_{qd}')
ylabel('Magnitude (dB)')
xticklabels({})
ylim([18 19])
grid on
hold off
subplot(2,1,2)
semilogx(freq_meas,angle(Z_qd_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,angle(Z_qd)./(2*pi)*360)
legend('Measured','Analytical','Location','southeast')
ylabel('Phase (deg)')
xlabel('Frequency (Hz)')
ylim([-5 5])
grid on
hold off

figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_qq_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,20.*log10(abs(Z_qq)))
legend('Measured','Analytical','Location','southeast')
title('Z_{qq}')
ylabel('Magnitude (dB)')
xticklabels({})
ylim([20 120])
grid on
hold off
subplot(2,1,2)
semilogx(freq_meas,angle(Z_qq_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,angle(Z_qq)./(2*pi)*360)
legend('Measured','Analytical','Location','southeast')
ylabel('Phase (deg)')
xlabel('Frequency (Hz)')
grid on
hold off
