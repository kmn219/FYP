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
                   
Y_sys_sym = simplify(inv(inv(Y_N)+inv(Y_A))); %inv(eye(3)+Y_N*inv(Y_A))*Y_N;
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
legend('Measured (voltage source) Z^{sys}_{1,1}','Theoretical Z^{sys}_{1,1}');
title('|Z^{sys}_{1,1}|')
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
grid on

figure
semilogx(V_freq_vals,vs_phase_11,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(I_freq_vals,is_phase_11,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f1/(2*pi),ph_11)
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
hold off
title('phase(Z^{sys}_{3,3})')
ylabel('Phase (°)')
xlabel('Frequency (Hz)')
grid on

%% Vector fitting

%close all
freq = logspace(-2,4,20);
vals_1 = freqresp(Z_sys_tf(3,3),2*pi.*freq);
clear vals
for i = 1:1:length(vals_1)
    vals(1,i) = vals_1(1,1,i);
end

opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=3;      %Include both D, E in fitting    
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %Do NOT skip identification of residues (C,D,E) 
opts.cmplx_ss=1;   %Create complex state space model

opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=1;       %Create magnitude plot for fitting of f(s) 
opts.logx=1;       %Use logarithmic abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis 
opts.errplot=1;    %Include deviation in magnitude plot
opts.phaseplot=0;  %Do not produce plot of phase angle (in addition to magnitiude)
opts.legend=1;     %Do include legends in plots

poles_sym = zeros(4,15);
poles_act = poles_sym;
poles_theo = poles_sym;

weight=ones(1,length(freq));

for i = [4]
    N=5+i;
    poles=-2*pi*logspace(log10(1),log10(1000),N);
%     bet=linspace(0,10,N/2);
%     poles=[];
%     for n=1:length(bet)
%       alf=-bet(n)*1e-2;
%       poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ];
%     end
    disp(['  N ' num2str(N)])
    Niter=10;
    for iter=1:Niter
        if iter==Niter, opts.skip_res=0; end
        %disp(['   Iter ' num2str(iter)])
        [SER,poles,rmserr,fit]=vectfit3(vals,2i*pi*freq,poles,weight,opts);
        rms(iter,1)=rmserr;
    end
    disp('Done Z_33.')
    poles_theo(i,1:length(poles)) = poles;
    SER_theo{i} = SER;
%     figure
%     title(['N: ' N])
%     scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf))
%     hold on
%     scatter(real(poles),imag(poles))
%     %ylim([-1 1])
%     %xlim([-2.5 0])
%     hold off
end

poles_theo(4,:)

%%


[res,p] = ss2pr(SER_theo{4}.A,SER_theo{4}.B,SER_theo{4}.C);
for i = 1:1:length(res)
    r(i,1) = res(1,1,i);
end
[b,a] = residue(r,p,[]);
v_fit_tf = tf(b,a);
clear error v_fit_mag theory_mag v_fit_m v_theory_m
v_fit_m = abs(freqresp(v_fit_tf,2*pi.*freq));
theory_m = abs(freqresp(Z_sys_tf(3,3),2*pi.*freq));
for i = 1:1:length(v_fit_m)
    v_fit_mag(1,i) = v_fit_m(1,1,i);
    theory_mag(1,i) = theory_m(1,1,i);
end
error = abs((v_fit_mag-theory_mag)./theory_mag).*100;

figure
semilogx(freq,error)
hold on
grid on
semilogx(freq,error,'LineStyle','none','Marker','x','MarkerSize',10,'MarkerEdgeColor',[0 0.4470 0.7410])
hold off
title('vectfit Approximation Error')
ylabel('Error (%)')
xlabel('Frequency (Hz)')
ylim([0 9e-13])

[m_b,p_b,f_b] = bode(v_fit_tf,{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});

clear magdb_b
for i = 1:1:length(m_b)
    magdb_b(1,i) = 20*log10(m_b(1,1,i));
end
for i = 1:1:length(p_b)
    ph_b(1,i) = p_b(1,1,i);
end

figure
semilogx(f_b./(2*pi),magdb_b,'Color','blue','LineWidth',1)
hold on
grid on
semilogx(freq,20*log10(abs(vals)),'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f3/(2*pi),magdb_33,'--','LineWidth',2,'Color','red')
title('vectfit Results with Nine Poles')
legend('vectfit model','Data used for vectfit','Theoretical model')


%%
f_sym = 10.^(magdb_33./20).*exp(1i*(ph_33./360.*2*pi));
weight=ones(1,length(f3));

for i = [1,2,3,4,5,10]
    N=5+i;
    poles=-2*pi*logspace(log10(1),log10(1000),N);
%     bet=linspace(0,10,N/2);
%     poles=[];
%     for n=1:length(bet)
%       alf=-bet(n)*1e-2;
%       poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ];
%     end
    disp(['  N ' num2str(N)])
    Niter=10;
    for iter=1:Niter
        if iter==Niter, opts.skip_res=0; end
        %disp(['   Iter ' num2str(iter)])
        [SER,poles,rmserr,fit]=vectfit3(f_sym,1i*f3,poles,weight,opts);
        rms(iter,1)=rmserr;
    end
    disp('Done Z_33.')
%     poles_sym(i,1:length(poles)) = poles;
%     SER_sym{i} = SER;
%     figure
%     scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf))
%     hold on
%     scatter(real(poles),imag(poles))
%     hold off
%     poles
end

clear f_act
f_act = 10.^(magnitude_db_33./20)'.*exp(1i*(phase_33./360.*2*pi)');
weight=ones(1,length(I_freq_vals));

for i = [1:1:5 10]
    N=5+i;
    poles=-2*pi*logspace(log10(1),log10(1000),N);
%     bet=linspace(0,10,N/2);
%     poles=[];
%     for n=1:length(bet)
%       alf=-bet(n)*1e-2;
%       poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ];
%     end
    disp('Vector fitting Z_33...')
    Niter=10;
    disp(['  N ' num2str(N)])
    for iter=1:Niter
        if iter==Niter, opts.skip_res=0; end
        [SER,poles,rmserr,fit]=vectfit3(f_act,2i*pi*I_freq_vals,poles,weight,opts);
        rms(iter,1)=rmserr;
    end
    disp('Done Z_33.')
    poles_act(i,1:length(poles)) = poles;
    SER_act{i} = SER;
    figure
    scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf))
    hold on
    scatter(real(poles),imag(poles))
    hold off
    poles
end

%%
figure
semilogx(I_freq_vals,20*log10(abs(f_act)))
hold on
semilogx(f3./(2*pi),magdb_33)

figure
plot(I_freq_vals,real(f_act))
hold on
plot(freq,real(vals))

%%
%weight=ones(1,length(f3));
weight=ones(1,length(I_freq_vals));

N=9;
poles=-2*pi*logspace(log10(1),log10(1000),N);
    disp(['  N ' num2str(N)])
    Niter=10;
%     for iter=1:Niter
%         if iter==Niter, opts.skip_res=0; end
%         %disp(['   Iter ' num2str(iter)])
%         [SER,poles,rmserr,fit]=vectfit3(f_sym,1i*f3,poles,weight,opts);
%         rms(iter,1)=rmserr;
%     end
    for iter=1:Niter
        if iter==Niter, opts.skip_res=0; end
        [SER,poles,rmserr,fit]=vectfit3(f_act,2i*pi*I_freq_vals,poles,weight,opts);
        rms(iter,1)=rmserr;
    end
    figure
    subplot(2,1,1)
    scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf))
    subplot(2,1,2)
    scatter(real(poles),imag(poles))
    poles

%%


[res,p] = ss2pr(SER_theo{4}.A,SER_theo{4}.B,SER_theo{4}.C);
for i = 1:1:length(res)
    r(i,1) = res(1,1,i);
end
[b,a] = residue(r,p,[]);
v_fit_tf = tf(b,a);
clear v_fit_mag theory_mag
v_fit_m = abs(freqresp(v_fit_tf,2*pi.*I_freq_vals));
theory_m = abs(freqresp(Z_sys_tf(3,3),2*pi.*I_freq_vals));
for i = 1:1:length(v_fit_m)
    v_fit_mag(1,i) = v_fit_m(1,1,i);
    theory_mag(1,i) = theory_m(1,1,i);
end
clear error
error = abs((v_fit_mag-theory_mag)./theory_mag).*100;

figure
semilogx(I_freq_vals,error)
hold on
grid on
semilogx(I_freq_vals,error,'LineStyle','none','Marker','x','MarkerSize',10,'MarkerEdgeColor',[0 0.4470 0.7410])
hold off
title('vectfit Approximation Error')
ylabel('Error (%)')
xlabel('Frequency (Hz)')

[m_b,p_b,f_b] = bode(v_fit_tf,{freq_lim(1,1)*2*pi,freq_lim(1,2)*2*pi});

clear magdb_b
for i = 1:1:length(m_b)
    magdb_b(1,i) = 20*log10(m_b(1,1,i));
end
for i = 1:1:length(p_b)
    ph_b(1,i) = p_b(1,1,i);
end

figure
semilogx(f_b./(2*pi),magdb_b,'Color','blue','LineWidth',1)
hold on
grid on
semilogx(I_freq_vals,magnitude_db_33,'LineStyle','none','Marker','x','MarkerSize',10)
semilogx(f3/(2*pi),magdb_33,'--','LineWidth',2,'Color','red')
title('vectfit Results with Nine Poles')
legend('vectfit model','Data used for vectfit','Theoretical model')

%%
figure
bode(tf(ss2tf(SER_act{3}.A,SER_act{3}.B,SER_act{3}.C,SER_act{3}.D)))

%%

figure
scatter(real(Z_sys_poles_tf),imag(Z_sys_poles_tf),'Marker','x','LineWidth',1,'SizeData',80)
hold on
scatter(real(poles_act(3,:)),imag((poles_act(3,:))),'Marker','x','LineWidth',1,'SizeData',80)
hold off
grid on
legend('tf poles','vfit poles (N=6)','vfit poles (N=7)','vfit poles (N=8)','vfit poles (N=9)')
xlabel('Real')
ylabel('Imaginary')
