%% Values definition
clear
clc
close all

freq_lim = [10^-2, 10^4];
tSettle = 1;

%% System Model
s = sym('s');

y12 = 1/(0.3*s+0.5);
y13 = 1/(0.2*s+0.6);
y23 = 1/(0.8*s+1.5);
y11 = 1/(5.8*s+1.2) + 7.9*s;
y22 = 1/(2*s+2.2) + 4*s;
y33 = 1/(5*s+5) + 6*s;

% y11 = tf([1],[5.8 1.2]) + tf([7.9 0],[1]);
% y12 = tf([1],[0.3 0.5]);
% y13 = tf([1],[0.2 0.6]);
% y23 = tf([1],[0.8 1.5]);
% y22 = tf([1],[2 2.2]) + tf([4 0],[1]);
% y33 = tf([1],[5 5]) + tf([6 0],[1]);


Y_nodal_sym = [y11+y12+y13, -y12, -y13
           -y12, y22+y12+y23, -y23
           -y13, -y23, y33+y23+y13];

Z_sys_sym = inv(Y_nodal_sym);
% Z_sys = minreal(inv(Y_nodal))

[N,D] = numden(Z_sys_sym(1,1));
Z_sys_poles_sym = vpasolve(D);

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

% figure
% scatter(real(zero(Y_nodal_tf)),imag(zero(Y_nodal_tf)))
% grid on

Z_sys_poles_tf = pole(Z_sys_tf(1,1));

figure
subplot(2,1,1)
scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf))
grid on
title('Poles of Z_{sys} (tf calculation)')
subplot(2,1,2)
scatter(real(Z_sys_poles_sym),imag(Z_sys_poles_sym))
grid on
title('Poles of Z_{sys} (sym calculation)')

Y_nodal_zeros = zero(Y_nodal_tf);

[m1,p1,f1] = bode(Z_sys_tf(1,1),{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});
[m2,p2,f2] = bode(Z_sys_tf(2,2),{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});
[m3,p3,f3] = bode(Z_sys_tf(3,3),{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});

%% Simulation and data gathering

V_amp = 20;
V_freq_vals = [0 logspace(log10(freq_lim(1,1)),log10(freq_lim(1,2)),19)];
for i = 1:1:length(V_freq_vals)
    tStep_temp = 1/(20*V_freq_vals(1,i));
    if tStep_temp > 1 | isinf(tStep_temp)
        tStep_vals(1,i) = 1;
    else
        tStep_vals(1,i) = tStep_temp;
    end
end
tStop_vals = [tSettle+1 tSettle+10./V_freq_vals(1,2:end)];

for i = 1:length(V_freq_vals)
    i
    tStop = tStop_vals(1,i);
    tStep = tStep_vals(1,i);
    V_freq = V_freq_vals(1,i);
    
    out = sim('basic_model_11.slx');
    len = length(out.ScopeData{1}.Values.Time(:,1));
    I_out_11(1:len,i) = out.ScopeData{2}.Values.Data(:,1);
    V_out_11(1:len,i) = out.ScopeData{1}.Values.Data(:,1);
    out = sim('basic_model_22.slx');
    I_out_22(1:len,i) = out.ScopeData{2}.Values.Data(:,1);
    V_out_22(1:len,i) = out.ScopeData{1}.Values.Data(:,1);
    out = sim('basic_model_33.slx');
    I_out_33(1:len,i) = out.ScopeData{2}.Values.Data(:,1);
    V_out_33(1:len,i) = out.ScopeData{1}.Values.Data(:,1);
    
    [A_I,~] = max(I_out_11(round(tSettle/tStep):len,i));
    [A_V,~] = max(V_out_11(round(tSettle/tStep):len,i));
    magnitude_db_11(i,1) = 20*log10(A_V/A_I);
    [A_I,~] = max(I_out_22(round(tSettle/tStep):len,i));
    [A_V,~] = max(V_out_22(round(tSettle/tStep):len,i));
    magnitude_db_22(i,1) = 20*log10(A_V/A_I);
    [A_I,~] = max(I_out_33(round(tSettle/tStep):len,i));
    [A_V,~] = max(V_out_33(round(tSettle/tStep):len,i));
    magnitude_db_33(i,1) = 20*log10(A_V/A_I);
end

%% Clear simulation values (large)
clear I_out_11 I_out_22 I_out_33 V_out_11 V_out_22 V_out_33 out

%% Impedance visualisation

for i = 1:1:length(m1)
    magdb_11(1,i) = 20*log10(m1(1,1,i));
end
for i = 1:1:length(m2)
    magdb_22(1,i) = 20*log10(m2(1,1,i));
end
for i = 1:1:length(m3)
    magdb_33(1,i) = 20*log10(m3(1,1,i));
end

figure
semilogx(V_freq_vals,magnitude_db_11,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(f1/(2*pi),magdb_11)
semilogx(V_freq_vals,magnitude_db_22,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f2/(2*pi),magdb_22)
semilogx(V_freq_vals,magnitude_db_33,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f3/(2*pi),magdb_33)
hold off
legend('Measured Z^{sys}_{1,1}','Theoretical Z^{sys}_{1,1}','Measured Z^{sys}_{2,2}','Theoretical Z^{sys}_{2,2}','Measured Z^{sys}_{3,3}','Theoretical Z^{sys}_{3,3}')
title('|Z^{sys}|')
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
grid on


%% Mode analysis

sz_Z_sys = size(Z_sys_sym);

% Implement eq 3.33
xi_eq = -trace(adjoint(Y_nodal_sym))/diff(det(Y_nodal_sym));
for i = 1:1:length(Z_sys_poles_sym)
    xi(i,1) = subs(xi_eq,s,Z_sys_poles_sym(i,1));
    res_Z_sys(1:sz_Z_sys(1,1),1:sz_Z_sys(1,2),i) = subs(adjoint(Y_nodal_sym)/diff(det(Y_nodal_sym)),s,Z_sys_poles_sym(i,1));
end

%% Residue calculation

residue(adjoint(Y_nodal_sym)/diff(det(Y_nodal_sym)))


