clear
close all

% Parameters
R_a = 15;           %Resistance (apparatus)
L_a = 27e-3;        %Inductance (apparatus)
R_g = 20;           %Resistance (grid)
L_g = 50e-3;        %Inductance (grid)
f_c = 50;           %Grid freq (Hz)
f_min = 10;         %Min measurement freq (Hz)
f_max = 10^5;       %Max measurement freq (Hz)
n = 20;             %No of measurement points
tStep_min = 5e-6;   %Sim step size
tSet = 50e-3;       %Settling time
amplitude = 100;    %Injection amplitude
n_cycles = 50;

doubleSidedMeasure = 1;     %1 = measure both sides

%% Calculate analytical results
f = sym('f');
Z_ana_a = [R_a+1i*2*pi*f*L_a, -2*pi*f_c*L_a
           2*pi*f_c*L_a,      R_a+1i*2*pi*f*L_a];
Z_ana_g = [R_g+1i*2*pi*f*L_g, -2*pi*f_c*L_g
           2*pi*f_c*L_g,      R_g+1i*2*pi*f*L_g];
% Z_ana = Z_ana_a;
Z_ana = inv(inv(Z_ana_a)+inv(Z_ana_g));

% s = tf('s');
% Z_sys = [R+s*L, -2*pi*f_c*L
%          2*pi*f_c*L, R+s*L];

f_max = min(f_max,1/tStep_min);

freq_vals = logspace(log10(f_min),log10(f_max),200);

% for i = 1:1:length(freq_vals)
%     Z_ana_a_val = eval(subs(Z_ana_a,f,freq_vals(1,i)));
%     Z_ana_g_val = eval(subs(Z_ana_g,f,freq_vals(1,i)));
%     Z_ana_val = inv(inv(Z_ana_a_val)+inv(Z_ana_g_val));
%     
%     Z_dd(1,i) = Z_ana_val(1,1);
%     Z_dq(1,i) = Z_ana_val(1,2);
%     Z_qd(1,i) = Z_ana_val(2,1);
%     Z_qq(1,i) = Z_ana_val(2,2);
% end

Z_dd = eval(subs(Z_ana(1,1),f,freq_vals));
Z_dq = eval(subs(Z_ana(1,2),f,freq_vals));
Z_qd = eval(subs(Z_ana(2,1),f,freq_vals));
Z_qq = eval(subs(Z_ana(2,2),f,freq_vals));

%% Run model and gather data

freq_meas = logspace(log10(f_min),log10(f_max),n);
for i = 1:1:length(freq_meas)
    i
    f_meas = freq_meas(1,i);
    t_meas = 1/f_meas;
    
    tStep = tStep_min;
    Fs = 1/tStep_min;
    tSettle = tStep_min*max(floor(tSet/t_meas),1);
    tStop = t_meas*n_cycles+tSettle;
    measBegin = floor(tSettle/tStep_min)+1;
    freq = [f_meas f_meas];
    phase = [1/2*pi 1/2*pi];
    
    % Run first simulation ================================================
    amp = [amplitude, 0];
    out = sim('dq0_double.slx');

    len = length(out.dq0{1}.Values.Data(measBegin:end,1));
    f = Fs*(0:(len/2))/len;

    % Apparatus side --------------------------------------------
    Y = fft(out.dq0{2}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd1_a_mag = abs(interp1(f,oneSide,f_meas));
    Vd1_a_ph = angle(interp1(f,oneSide,f_meas));
    Vd1_a = Vd1_a_mag*exp(1i*0);

    Y = fft(out.dq0{2}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq1_a_mag = abs(interp1(f,oneSide,f_meas));
    Vq1_a_ph = angle(interp1(f,oneSide,f_meas))-Vd1_a_ph;
    Vq1_a = Vq1_a_mag*exp(1i*Vq1_a_ph);

    Y = fft(out.dq0{4}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id1_a_mag = abs(interp1(f,oneSide,f_meas));
    Id1_a_ph = angle(interp1(f,oneSide,f_meas))-Vd1_a_ph;
    Id1_a = Id1_a_mag*exp(1i*Id1_a_ph);

    Y = fft(out.dq0{4}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq1_a_mag = abs(interp1(f,oneSide,f_meas));
    Iq1_a_ph = angle(interp1(f,oneSide,f_meas))-Vd1_a_ph;
    Iq1_a = Iq1_a_mag*exp(1i*Iq1_a_ph);

    % Grid side -------------------------------------------------
    Y = fft(out.dq0{3}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd1_g_mag = abs(interp1(f,oneSide,f_meas));
    Vd1_g_ph = angle(interp1(f,oneSide,f_meas));
    Vd1_g = Vd1_g_mag*exp(1i*0);

    Y = fft(out.dq0{3}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq1_g_mag = abs(interp1(f,oneSide,f_meas));
    Vq1_g_ph = angle(interp1(f,oneSide,f_meas))-Vd1_g_ph;
    Vq1_g = Vq1_g_mag*exp(1i*Vq1_g_ph);

    Y = fft(out.dq0{5}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id1_g_mag = abs(interp1(f,oneSide,f_meas));
    Id1_g_ph = angle(interp1(f,oneSide,f_meas))-Vd1_g_ph;
    Id1_g = Id1_g_mag*exp(1i*Id1_g_ph);

    Y = fft(out.dq0{5}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq1_g_mag = abs(interp1(f,oneSide,f_meas));
    Iq1_g_ph = angle(interp1(f,oneSide,f_meas))-Vd1_g_ph;
    Iq1_g = Iq1_g_mag*exp(1i*Iq1_g_ph);
    
    % Run second simulation ===============================================
    clear out dq0
    amp = [0, amplitude];
    out = sim('dq0_double.slx');

    % Apparatus side --------------------------------------------
    Y = fft(out.dq0{2}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq2_a_mag = abs(interp1(f,oneSide,f_meas));
    Vq2_a_ph = angle(interp1(f,oneSide,f_meas));
    Vq2_a = Vq2_a_mag*exp(1i*0);

    Y = fft(out.dq0{2}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd2_a_mag = abs(interp1(f,oneSide,f_meas));
    Vd2_a_ph = angle(interp1(f,oneSide,f_meas))-Vq2_a_ph;
    Vd2_a = Vd2_a_mag*exp(1i*Vd2_a_ph);

    Y = fft(out.dq0{4}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id2_a_mag = abs(interp1(f,oneSide,f_meas));
    Id2_a_ph = angle(interp1(f,oneSide,f_meas))-Vq2_a_ph;
    Id2_a = Id2_a_mag*exp(1i*Id2_a_ph);

    Y = fft(out.dq0{4}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Iq2_a_mag = abs(interp1(f,oneSide,f_meas));
    Iq2_a_ph = angle(interp1(f,oneSide,f_meas))-Vq2_a_ph;
    Iq2_a = Iq2_a_mag*exp(1i*Iq2_a_ph);

    % Grid side -------------------------------------------------
    Y = fft(out.dq0{3}.Values.Data(measBegin:end,2));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vq2_g_mag = abs(interp1(f,oneSide,f_meas));
    Vq2_g_ph = angle(interp1(f,oneSide,f_meas));
    Vq2_g = Vq2_g_mag*exp(1i*0);

    Y = fft(out.dq0{3}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Vd2_g_mag = abs(interp1(f,oneSide,f_meas));
    Vd2_g_ph = angle(interp1(f,oneSide,f_meas))-Vq2_g_ph;
    Vd2_g = Vd2_g_mag*exp(1i*Vd2_g_ph);

    Y = fft(out.dq0{5}.Values.Data(measBegin:end,1));
    oneSide = 2*Y(1:floor(len/2)+1)/len;
    Id2_g_mag = abs(interp1(f,oneSide,f_meas));
    Id2_g_ph = angle(interp1(f,oneSide,f_meas))-Vq2_g_ph;
    Id2_g = Id2_g_mag*exp(1i*Id2_g_ph);

    Y = fft(out.dq0{5}.Values.Data(measBegin:end,2));
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
        Z_meas = inv(inv(Z_a)+inv(Z_g))
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

figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_dd_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,20.*log10(abs(Z_dd)))
title('Z_{dd}')
ylabel('Magnitude (dB)')
legend('Measured','Analytical','Location','southeast')
xticklabels({})
ylim([0 100])
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

figure
subplot(2,1,1)
semilogx(freq_meas,20.*log10(abs(Z_dq_meas)),'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(freq_vals,20.*log10(abs(Z_dq)))
legend('Measured','Analytical','Location','southeast')
title('Z_{dq}')
ylabel('Magnitude (dB)')
xticklabels({})
ylim([14 16])
grid on
hold off
subplot(2,1,2)
angles = angle(Z_dq_meas)./(2*pi)*360;
for i = 1:1:length(angles)
    if angles(1,i) < 0
        angles(1,i) = -angles(1,i);
    end
end
semilogx(freq_meas,angles,'LineStyle','none','Marker','x','MarkerSize',10)
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
ylim([14 16])
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
ylim([0 100])
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
