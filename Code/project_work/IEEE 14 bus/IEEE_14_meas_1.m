clear
close all

% Parameters
f_c = 50;           %Grid freq (Hz)
f_min = 10^1;         %Min measurement freq (Hz)
f_max = 40;       %Max measurement freq (Hz)
n = 20;             %No of measurement points
Ts = 4e-5;   %Sim min step size
amplitude = 5/100;    %Injection amplitude
tSet = 10;
n_cycles = 20;
doubleSidedMeasure = 1;

load('parameters.mat')
load('initial_conds.mat')

%% Simulate

freq_meas = logspace(log10(f_min),log10(f_max),n);
for i = 1:1:length(freq_meas)
    i
    f_meas = freq_meas(1,i);
    t_meas = 1/f_meas;
    
    Fs = 1/Ts;
    tSettle = Ts*max(floor(tSet/t_meas),1);
    tStop = 100+t_meas*n_cycles+tSettle;
    measBegin = floor(tSettle/Ts)+1;
    freq = [f_meas f_meas];
    phase = [1/2*pi 1/2*pi];
    
    % Run first simulation ================================================
    amp = [amplitude, 0];
    out = sim('IEEE14.slx');

    len = length(ScopeData2{1}.Values.Data(measBegin:end,1));
    f = Fs*(0:(len/2))/len;

    % Apparatus side --------------------------------------------
    Y = fft(ScopeData2{2}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd1_a_mag = abs(interp1(f,oneSide,f_meas));
    Vd1_a_ph = angle(interp1(f,oneSide,f_meas));
    Vd1_a = Vd1_a_mag*exp(1i*0);

    Y = fft(ScopeData2{2}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq1_a_mag = abs(interp1(f,oneSide,f_meas));
    Vq1_a_ph = angle(interp1(f,oneSide,f_meas))-Vd1_a_ph;
    Vq1_a = Vq1_a_mag*exp(1i*Vq1_a_ph);

    Y = fft(ScopeData2{4}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id1_a_mag = abs(interp1(f,oneSide,f_meas));
    Id1_a_ph = angle(interp1(f,oneSide,f_meas))-Vd1_a_ph;
    Id1_a = Id1_a_mag*exp(1i*Id1_a_ph);

    Y = fft(ScopeData2{4}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq1_a_mag = abs(interp1(f,oneSide,f_meas));
    Iq1_a_ph = angle(interp1(f,oneSide,f_meas))-Vd1_a_ph;
    Iq1_a = Iq1_a_mag*exp(1i*Iq1_a_ph);

    % Grid side -------------------------------------------------
    Y = fft(ScopeData2{3}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd1_g_mag = abs(interp1(f,oneSide,f_meas));
    Vd1_g_ph = angle(interp1(f,oneSide,f_meas));
    Vd1_g = Vd1_g_mag*exp(1i*0);

    Y = fft(ScopeData2{3}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq1_g_mag = abs(interp1(f,oneSide,f_meas));
    Vq1_g_ph = angle(interp1(f,oneSide,f_meas))-Vd1_g_ph;
    Vq1_g = Vq1_g_mag*exp(1i*Vq1_g_ph);

    Y = fft(ScopeData2{5}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id1_g_mag = abs(interp1(f,oneSide,f_meas));
    Id1_g_ph = angle(interp1(f,oneSide,f_meas))-Vd1_g_ph;
    Id1_g = Id1_g_mag*exp(1i*Id1_g_ph);

    Y = fft(ScopeData2{5}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq1_g_mag = abs(interp1(f,oneSide,f_meas));
    Iq1_g_ph = angle(interp1(f,oneSide,f_meas))-Vd1_g_ph;
    Iq1_g = Iq1_g_mag*exp(1i*Iq1_g_ph);
    
    % Run second simulation ===============================================
    clear out dq0
    amp = [0, amplitude];
    out = sim('IEEE14.slx');

    % Apparatus side --------------------------------------------
    Y = fft(ScopeData2{2}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq2_a_mag = abs(interp1(f,oneSide,f_meas));
    Vq2_a_ph = angle(interp1(f,oneSide,f_meas));
    Vq2_a = Vq2_a_mag*exp(1i*0);

    Y = fft(ScopeData2{2}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd2_a_mag = abs(interp1(f,oneSide,f_meas));
    Vd2_a_ph = angle(interp1(f,oneSide,f_meas))-Vq2_a_ph;
    Vd2_a = Vd2_a_mag*exp(1i*Vd2_a_ph);

    Y = fft(ScopeData2{4}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id2_a_mag = abs(interp1(f,oneSide,f_meas));
    Id2_a_ph = angle(interp1(f,oneSide,f_meas))-Vq2_a_ph;
    Id2_a = Id2_a_mag*exp(1i*Id2_a_ph);

    Y = fft(ScopeData2{4}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq2_a_mag = abs(interp1(f,oneSide,f_meas));
    Iq2_a_ph = angle(interp1(f,oneSide,f_meas))-Vq2_a_ph;
    Iq2_a = Iq2_a_mag*exp(1i*Iq2_a_ph);

    % Grid side -------------------------------------------------
    Y = fft(ScopeData2{3}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq2_g_mag = abs(interp1(f,oneSide,f_meas));
    Vq2_g_ph = angle(interp1(f,oneSide,f_meas));
    Vq2_g = Vq2_g_mag*exp(1i*0);

    Y = fft(ScopeData2{3}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd2_g_mag = abs(interp1(f,oneSide,f_meas));
    Vd2_g_ph = angle(interp1(f,oneSide,f_meas))-Vq2_g_ph;
    Vd2_g = Vd2_g_mag*exp(1i*Vd2_g_ph);

    Y = fft(ScopeData2{5}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id2_g_mag = abs(interp1(f,oneSide,f_meas));
    Id2_g_ph = angle(interp1(f,oneSide,f_meas))-Vq2_g_ph;
    Id2_g = Id2_g_mag*exp(1i*Id2_g_ph);

    Y = fft(ScopeData2{5}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq2_g_mag = abs(interp1(f,oneSide,f_meas));
    Iq2_g_ph = angle(interp1(f,oneSide,f_meas))-Vq2_g_ph;
    Iq2_g = Iq2_g_mag*exp(1i*Iq2_g_ph);

    % Calculate impedance --------------------------------------
    vs_a = [Vd1_a 0
            0     Vq2_a];
    is_a = [Id1_a Id2_a
            Iq1_a Iq2_a];
    Z_a = inv(is_a)*vs_a;

    if doubleSidedMeasure == 1
        vs_g = [Vd1_g 0
                0     Vq2_g];
        is_g = [Id1_g Id2_g
                Iq1_g Iq2_g];
        Z_g = inv(is_g)*vs_g;
        Z_meas = inv(inv(Z_a)+inv(Z_g));
    else
        Z_meas = Z_a;
    end
    Z_dd_meas(1,i) = Z_meas(1,1);
    Z_dq_meas(1,i) = Z_meas(1,2);
    Z_qd_meas(1,i) = Z_meas(2,1);
    Z_qq_meas(1,i) = Z_meas(2,2);
end

%% Plot data
close all
clear

load('measurement1.mat')
load('measurement2.mat')
load('Z1212.mat')

for i = 1:1:length(Z1212_mag)
    Z1212_mag_new(1,i) = Z1212_mag(1,1,i);
    if Z1212_phase(1,1,i) < -360
        Z1212_phase_new(1,i) = Z1212_phase(1,1,i)+4*360;
    else
        Z1212_phase_new(1,i) = Z1212_phase(1,1,i);
    end
end



figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_dd_meas)),'LineStyle','none','Marker','x','MarkerSize',10,'MarkerEdgeColor',[0 0.4470 0.7410])
hold on
semilogx(Z1212_wout/(2*pi),20.*log10(Z1212_mag_new),'Color',[0.8500 0.3250 0.0980])
semilogx(freq_meas1,20.*log10(abs(Z_dd_meas1)),'LineStyle','none','Marker','x','MarkerSize',10,'MarkerEdgeColor',[0 0.4470 0.7410])
xlim([10^0 10^3])
title('Z^{sys}_{12,12} in dd axis')
ylabel('Magnitude (dB)')
legend('Measurement','Analytical','Location','southeast')
xticklabels({})
ylim([-25 15])
grid on
hold off
subplot(2,1,2)
semilogx(freq_meas,angle(Z_dd_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10,'MarkerEdgeColor',[0 0.4470 0.7410])
hold on
semilogx(Z1212_wout/(2*pi),Z1212_phase_new,'Color',[0.8500 0.3250 0.0980])
semilogx(freq_meas1,angle(Z_dd_meas1)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10,'MarkerEdgeColor',[0 0.4470 0.7410])
xlim([10^0 10^3])
legend('Measurement','Analytical','Location','southeast')
ylabel('Phase (deg)')
xlabel('Frequency (Hz)')
grid on
hold off

% Z_dq_meas(1,18) = conj(Z_dq_meas(1,18));
% Z_dq_meas(1,19) = conj(Z_dq_meas(1,19));
% 
% figure
% subplot(2,1,1)
% semilogx(freq_meas,20.*log10(abs(Z_dq_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
% hold on
% semilogx(freq_vals,20.*log10(abs(Z_dq)))
% legend('Measured','Analytical','Location','southeast')
% title('Z_{dq}')
% ylabel('Magnitude (dB)')
% xticklabels({})
% ylim([18 19])
% grid on
% hold off
% subplot(2,1,2)
% semilogx(freq_meas,angle(Z_dq_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
% hold on
% semilogx(freq_vals,angle(Z_dq)./(2*pi)*360)
% legend('Measured','Analytical','Location','southeast')
% ylabel('Phase (deg)')
% xlabel('Frequency (Hz)')
% ylim([175 185])
% grid on
% hold off
% 
% figure
% subplot(2,1,1)
% semilogx(freq_meas,20.*log10(abs(Z_qd_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
% hold on
% semilogx(freq_vals,20.*log10(abs(Z_qd)))
% legend('Measured','Analytical','Location','southeast')
% title('Z_{qd}')
% ylabel('Magnitude (dB)')
% xticklabels({})
% ylim([18 19])
% grid on
% hold off
% subplot(2,1,2)
% semilogx(freq_meas,angle(Z_qd_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
% hold on
% semilogx(freq_vals,angle(Z_qd)./(2*pi)*360)
% legend('Measured','Analytical','Location','southeast')
% ylabel('Phase (deg)')
% xlabel('Frequency (Hz)')
% ylim([-5 5])
% grid on
% hold off
% 
% figure
% subplot(2,1,1)
% semilogx(freq_meas,20.*log10(abs(Z_qq_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
% hold on
% semilogx(freq_vals,20.*log10(abs(Z_qq)))
% legend('Measured','Analytical','Location','southeast')
% title('Z_{qq}')
% ylabel('Magnitude (dB)')
% xticklabels({})
% ylim([20 120])
% grid on
% hold off
% subplot(2,1,2)
% semilogx(freq_meas,angle(Z_qq_meas)./(2*pi)*360,'LineStyle','none','Marker','x','MarkerSize',10)
% hold on
% semilogx(freq_vals,angle(Z_qq)./(2*pi)*360)
% legend('Measured','Analytical','Location','southeast')
% ylabel('Phase (deg)')
% xlabel('Frequency (Hz)')
% grid on
% hold off