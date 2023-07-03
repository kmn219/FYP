%% Run simulation
clear
close all

V_DC_vals = [1 zeros(1,9)];
V_AC_amp = 1;
V_AC_freq_vals = [0 logspace(0,3,9)];
tStop_vals = [7 4+10./V_AC_freq_vals(1,2:end)];
tStep_vals = [1e-3 1./(20.*V_AC_freq_vals(1,2:end))];
R = 10;
L = 1;

for i = 1:length(V_DC_vals)
    i
    tStop = tStop_vals(1,i);
    tStep = tStep_vals(1,i);
    V_DC = V_DC_vals(1,i);
    V_AC_freq = V_AC_freq_vals(1,i);
    t_temp = 0:tStep:tStop;
    len(1,i) = length(t_temp);
    t(1:len(1,i),i) = t_temp';
    V(1:len(1,i),i) = (V_DC + V_AC_amp.*sin(V_AC_freq.*t_temp))';
    out = sim('simple_model.slx');
%     time = out.scope_Iout{1}.Values.Time;
    I_out(1:len(1,i),i) = out.scope_Iout{1}.Values.Data(:,1);
    data(1:len(1,i),i) = I_out(1:len(1,i),i)./V(1:len(1,i),i);
end
%% Plot signals

figure
hold on
for i = 1:length(V_DC_vals)
    plot(t(1:len(1,i),i),I_out(1:len(1,i),i));
end
hold off
xlim([3 4]);
ylabel('I_{out} (A)')
grid on

figure
hold on
for i = 1:length(V_DC_vals)
    plot(t(1:len(1,i),i),V(1:len(1,i),i));
end
xlim([3 4]);
ylabel('V_{in} (V)')
grid on

%% Get data points and plot 

close all
t_periods = 1./V_AC_freq_vals;
for i = 1:length(V_DC_vals)
    tStep = tStep_vals(1,i); 
    [A,I] = max(I_out(round(4/tStep):len(1,i),i));
    magnitude_db(i,1) = 20*log10(A);
    index = rem(I,20);
    phase = 0;
end

[mag,phase,wout] = bode(tf([1],[L R]));
magdb_temp = 20*log10(mag);

for i = 1:1:length(magdb_temp)
    magdb(i,1) = magdb_temp(1,1,i);
end

figure
semilogx(V_AC_freq_vals*2*pi,magnitude_db,'LineStyle','none','Marker','x','MarkerSize',10)
hold on
semilogx(wout,magdb)
hold off

figure
plot(t(:,9),I_out(:,9),t(:,9),I_out)
grid on