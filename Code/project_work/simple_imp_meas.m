%% Define Values
clear
clc
close all

freq_lim = [10^-2, 10^2];
tSettle = 0.25;
V_amp = 20;
I_amp = 20;
n = 10;

% Options:
plotPoles = 0;

%% Define System Model
s = sym('s');

y12 = 1/(0.3*s+0.5);
y13 = 1/(0.2*s+0.6);
y23 = 1/(0.8*s+1.5);
y1 = 1/(5.8*s+1.2) + 7.9*s;
y2 = 1/(2*s+2.2) + 4*s;
y3 = 1/(5*s+5) + 6*s;

Y_N = [0    -y12 -y13
       -y12  0   -y23
       -y13 -y23  0  ];
   
Y_A = [y1+y12+y13 0          0
       0          y2+y12+y23 0
       0          0          y3+y13+y23];

Y_A = simplify(Y_A);
Y_N = simplify(Y_N);

Y_nodal_sym = [y1+y12+y13, -y12, -y13
              -y12, y2+y12+y23, -y23
              -y13, -y23, y3+y23+y13];
                   
Y_sys_sym = simplify(inv(inv(Y_N)+inv(Y_A)));
Z_sys_sym = inv(Y_nodal_sym);

%% Pole/zero calculation
close all

sz = size(Z_sys_sym);

for i = 1:1:sz(1,1)
    for j = 1:1:sz(1,2)
        [N,D] = numden(Z_sys_sym(i,j));
        N = sym2poly(N);
        D = sym2poly(D);
        Z_sys_tf(i,j) = tf(N,D);
    end
end

for i = 1:1:sz(1,1)
    for j = 1:1:sz(1,2)
        [N,D] = numden(Y_nodal_sym(i,j));
        N = sym2poly(N);
        D = sym2poly(D);
        Y_nodal_tf(i,j) = tf(N,D);
    end
end

for i = 1:1:sz(1,1)
    for j = 1:1:sz(1,2)
        [N,D] = numden(Y_sys_sym(i,j));
        N = sym2poly(N);
        D = sym2poly(D);
        Y_sys_tf(i,j) = tf(N,D);
    end
end

[N,D] = numden(simplify(Y_sys_sym(1,1)));
Y_sys_poles_sym = vpasolve(D);

Z_sys_poles_tf = pole(Z_sys_tf(1,1));
Y_sys_poles_tf = pole(Y_sys_tf(1,1));

if plotPoles == 1
    figure
    subplot(2,2,1)
    scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf))
    grid on
    title('Poles of Z^{sys} (tf calculation)')
    subplot(2,2,3)
    scatter(real(Z_sys_poles_sym),imag(Z_sys_poles_sym))
    grid on
    title('Poles of Z^{sys} (sym calculation)')
    subplot(2,2,2)
    scatter(real(Y_sys_poles_sym),imag(Y_sys_poles_sym))
    grid on
    title('Poles of Y^{sys} (sym calculation)')
    subplot(2,2,4)
    scatter(real(Y_sys_poles_tf),imag(Y_sys_poles_tf))
    grid on
    title('Poles of Y^{sys} (tf calculation)')
end

[m1,p1,f1] = bode(Z_sys_tf(1,1),{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});
[m2,p2,f2] = bode(Z_sys_tf(2,2),{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});
[m3,p3,f3] = bode(Z_sys_tf(3,3),{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});

%% Voltage source simulation

V_freq_vals = [logspace(log10(freq_lim(1,1)),log10(freq_lim(1,2)),n)];
for i = 1:1:length(V_freq_vals)
    tStep_temp = 1/(20*V_freq_vals(1,i));
    if tStep_temp > 1 | isinf(tStep_temp)
        tStep_vals(1,i) = 1;
    else
        tStep_vals(1,i) = tStep_temp;
    end
end
tStop_vals = [tSettle+10./V_freq_vals];

for i = 1:length(V_freq_vals)
    i
    tStop = tStop_vals(1,i);
    tStep = tStep_vals(1,i);
    V_freq = V_freq_vals(1,i);
    
    out = sim('basic_model_11_vscm.slx');
    V_g_11 = out.V_g_mag(end,1)*exp(1i*out.V_g_phase(end,1)/360*2*pi);
    I_g_11 = out.I_g_mag(end,1)*exp(1i*out.I_g_phase(end,1)/360*2*pi);
    V_a_11 = V_g_11 - V_amp;
    I_a_11 = out.I_g_mag(end,1)*exp(1i*(out.I_g_phase(end,1)-180)/360*2*pi);
    out = sim('basic_model_22_vscm.slx');
    V_g_22 = out.V_g_mag(end,1)*exp(1i*out.V_g_phase(end,1)/360*2*pi);
    I_g_22 = out.I_g_mag(end,1)*exp(1i*out.I_g_phase(end,1)/360*2*pi);
    V_a_22 = V_g_22 - V_amp;
    I_a_22 = out.I_g_mag(end,1)*exp(1i*(out.I_g_phase(end,1)-180)/360*2*pi);
    out = sim('basic_model_33_vscm.slx');
    V_g_33 = out.V_g_mag(end,1)*exp(1i*out.V_g_phase(end,1)/360*2*pi);
    I_g_33 = out.I_g_mag(end,1)*exp(1i*out.I_g_phase(end,1)/360*2*pi);
%     V_a_33 = out.V_a_mag(end,1)*exp(1i*out.V_a_phase(end,1)/360*2*pi);
    V_a_33 = V_g_33 - V_amp;
%     I_a_33 = out.I_a_mag(end,1)*exp(1i*out.I_a_phase(end,1)/360*2*pi);
    I_a_33 = out.I_g_mag(end,1)*exp(1i*(out.I_g_phase(end,1)-180)/360*2*pi);
    
    Z11_meas = 1/(I_g_11/V_g_11+I_a_11/V_a_11);
    vs_magnitude_db_11(1,i) = 20.*log10(abs(Z11_meas));
    vs_phase_11(1,i) = angle(Z11_meas)/(2*pi)*360;
    Z22_meas = 1/(I_g_22/V_g_22+I_a_22/V_a_22);
    vs_magnitude_db_22(1,i) = 20.*log10(abs(Z22_meas));
    vs_phase_22(1,i) = angle(Z22_meas)/(2*pi)*360;
    Z33_meas = 1/(I_g_33/V_g_33+I_a_33/V_a_33);
%     Z33_pred = eval(subs(Z_sys_sym(3,3),s,2*pi*1i*V_freq))
    vs_magnitude_db_33(1,i) = 20.*log10(abs(Z33_meas));
    vs_phase_33(1,i) = angle(Z33_meas)/(2*pi)*360;
end

%% Current source simulation

freq_diff = (log10(freq_lim(1,2)) - log10(freq_lim(1,1)))/(n-1)/2;
I_freq_vals = [logspace(log10(freq_lim(1,1))+freq_diff,log10(freq_lim(1,2))-freq_diff,n-1)];
for i = 1:1:length(I_freq_vals)
    tStep_temp = 1/(20*I_freq_vals(1,i));
    if tStep_temp > 1 | isinf(tStep_temp)
        tStep_vals(1,i) = 1;
    else
        tStep_vals(1,i) = tStep_temp;
    end
end
tStop_vals = [tSettle+10./I_freq_vals];

for i = 1:length(I_freq_vals)
    i
    tStop = tStop_vals(1,i);
    tStep = tStep_vals(1,i);
    I_freq = I_freq_vals(1,i);
    
    out = sim('basic_model_11_csvm.slx');
    V_amp_11 = out.V_out_mag(end,1);
    V_phase_11 = out.V_out_phase(end,1);
    out = sim('basic_model_22_csvm.slx');
    V_amp_22 = out.V_out_mag(end,1);
    V_phase_22 = out.V_out_phase(end,1);
    out = sim('basic_model_33_csvm.slx');
    V_amp_33 = out.V_out_mag(end,1);
    V_phase_33 = out.V_out_phase(end,1);
    
    is_magnitude_db_11(i,1) = 20*log10(V_amp_11/I_amp);
    is_phase_11(i,1) = V_phase_11;
    is_magnitude_db_22(i,1) = 20*log10(V_amp_22/I_amp);
    is_phase_22(i,1) = V_phase_22;
    is_magnitude_db_33(i,1) = 20*log10(V_amp_33/I_amp);
    is_phase_33(i,1) = V_phase_33;
end

%% Clear simulation values (large)
clear out

%% Impedance visualisation 

for i = 1:1:length(m1)
    magdb_11(1,i) = 20.*log10(m1(1,1,i));
end
for i = 1:1:length(m2)
    magdb_22(1,i) = 20.*log10(m2(1,1,i));
end
for i = 1:1:length(m3)
    magdb_33(1,i) = 20.*log10(m3(1,1,i));
end
for i = 1:1:length(p1)
    ph_11(1,i) = p1(1,1,i);
end
for i = 1:1:length(p2)
    ph_22(1,i) = p2(1,1,i);
end
for i = 1:1:length(p3)
    ph_33(1,i) = p3(1,1,i);
end

figure
semilogx(V_freq_vals,vs_magnitude_db_11,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_magnitude_db_11,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f1/(2*pi),magdb_11)
hold off
legend('Measured (voltage source)','Measured (current source)','Theoretical');
title('|Z^{sys}_{1,1}|')
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
grid on

figure
semilogx(V_freq_vals,vs_phase_11,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_phase_11,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f1/(2*pi),ph_11)
legend('Measured (voltage source)','Measured (current source)','Theoretical');
hold off
title('phase(Z^{sys}_{1,1})')
ylabel('Phase (°)')
xlabel('Frequency (Hz)')
grid on

figure
semilogx(V_freq_vals,vs_magnitude_db_22,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_magnitude_db_22,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f2/(2*pi),magdb_22)
legend('Measured (voltage source)','Measured (current source)','Theoretical');
hold off
legend('Measured Z^{sys}_{2,2}','Theoretical Z^{sys}_{2,2}');
title('|Z^{sys}_{2,2}|')
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
grid on

figure
semilogx(V_freq_vals,vs_phase_22,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_phase_22,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f2/(2*pi),ph_22)
legend('Measured (voltage source)','Measured (current source)','Theoretical');
hold off
title('phase(Z^{sys}_{2,2})')
ylabel('Phase (°)')
xlabel('Frequency (Hz)')
grid on

figure
semilogx(V_freq_vals,vs_magnitude_db_33,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_magnitude_db_33,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f3/(2*pi),magdb_33)
legend('Measured (voltage source)','Measured (current source)','Theoretical');
hold off
legend('Measured Z^{sys}_{3,3}','Theoretical Z^{sys}_{3,3}');
title('|Z^{sys}_{3,3}|')
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
grid on

figure
semilogx(V_freq_vals,vs_phase_33,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_phase_33,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f3/(2*pi),ph_33)
legend('Measured (voltage source)','Measured (current source)','Theoretical');
hold off
title('phase(Z^{sys}_{3,3})')
ylabel('Phase (°)')
xlabel('Frequency (Hz)')
grid on