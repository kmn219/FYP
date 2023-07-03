%% data load
clear;
load('Zdetuned1.mat');
load('GROUP3. 15V-5V-7V-10V.mat')
Length=length(Bode_0)/4;
f_s=Bode_0(2,1:2:Length-1);
Ns=length(f_s);
%% extract frequency response
%Ns=500;
%f_s=linspace(1,1000,Ns); %Vector holding the frequency samples, [rad/sec]. 1Hz~1000Hz
%f_s=(1:1:1000);
w_s=2*pi*f_s;
Z12dd = freqresp( Zdetuned(23,23),f_s,'Hz'); Z12dd=reshape(Z12dd,[1,Ns]);
Z12dq = freqresp( Zdetuned(23,24),f_s,'Hz'); Z12dq=reshape(Z12dq,[1,Ns]);
Z12qd = freqresp( Zdetuned(24,23),f_s,'Hz'); Z12qd=reshape(Z12qd,[1,Ns]);
Z12qq = freqresp( Zdetuned(24,24),f_s,'Hz'); Z12qq=reshape(Z12qq,[1,Ns]);
Z13dd = freqresp( Zdetuned(25,25),f_s,'Hz'); Z13dd=reshape(Z13dd,[1,Ns]);

Zsys=[Z12dd;Z13dd];
%Zsys=[Z12dd;Z12dq;Z12qd;Z12qq];

%% vecotr fitting
%Initial poles for Vector Fitting:
N=20; %order of approximation

%Initial poles: Complex conjugate pairs, linearly spaced:
bet=linspace(w_s(1),w_s(end),N/2);
poles=[];
for k=1:length(bet)
alf=-bet(k)*1e-2;
poles=[poles (alf-1i*bet(k)) (alf+1i*bet(k)) ];
end
%poles=-2*pi*logspace(0,4,N); 
weight=ones(1,Ns); %No weight
%weight=zeros(1,Ns); %Strong inverse weight
%for k=1:Ns
%weight(1,k)=1/sqrt(norm(f_mea(:,k)));
%end
opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=1;      %Include both D, E in fitting    
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %Do NOT skip identification of residues (C,D,E) 
opts.cmplx_ss=1;   %Create complex state space model
opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=1;       %Create magnitude plot for fitting of f(s) 
opts.logx=1;       %Use logarithmic abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis 
opts.errplot=1;    %Include deviation in magnitude plot
opts.phaseplot=1;  %Also produce plot of phase angle (in addition to magnitiude)
opts.legend=1;     %include legends in plots
%opts.cmplx_ss=1;



%% Forming (weighted) column sum:
g=0;
Nc=1;
for n=1:Nc
  %g=g+f(n,:); %unweighted sum     
  g=g+Zsys(n,:);%/norm(Zsys(n,:));
  %g=g+f(n,:)/sqrt(norm(f(n,:)));     
end
%weight_g=1./abs(g);
weight_g=ones(1,length(g));

disp('****Calculating improved initial poles by fitting column sum ...')
Niter1=3;
for iter=1:Niter1
   disp(['   Iter ' num2str(iter)])
   if iter==Niter1,opts.skip_res=0; end
   [SER,poles,rmserr,fit]=vectfit3(g,1i*w_s,poles,weight_g,opts);  
end

%% fitting
disp('vector fitting...')
%opts.skip_res=1;
weight=1./sqrt(abs(Zsys));
for j=1:3
    [SER,poles,rmserr,fit]=vectfit3(Zsys,1i*w_s,poles,weight,opts); 
end

AA=[];
BB=[];
CC=[];
tell=0;
for row=1:2
    AA=blkdiag(AA,SER.A);
    BB=blkdiag(BB,SER.B);
    for col=2:2
        tell=tell+1;
        CC(row,(col-1)*N+1:col*N)=SER.C(tell,:); 
        CC(col,(row-1)*N+1:row*N)=SER.C(tell,:);   
    end
end
%SER=col2full(SER);
%SER=tri2full(SER);
[Zsys_res,Zsys_pole] = ss2pr(AA,BB,CC);
Zsys_pole_Hz=Zsys_pole/2/pi;
disp('Done.')
Zsys_ss=ss(SER.A, SER.B, SER.C, SER.D);
Zsys12_fit=tf(Zsys_ss(1));
[num,den] = tfdata(Zsys12_fit); num=num{1}; den=den{1};
[r,p,k] = residue(num,den); p=p/2/pi;
figure(11)
clf;
bode(Zsys_ss(1)); hold on;
bode(Zsys_ss(2)); hold on;
%clf
%bode(SER.A, SER.B, SER.C, SER.D)
%Res_dd=reshape(Zsys_res(1,1,:),[length(Zsys_res),1]);


%% plot to verify
Ydd_mag_db = 20*log10(abs(Z12dd));
Ydd_Phase = angle(Z12dd)/2/pi*360; 
Ydd_mag_db2 = 20*log10(abs(Z13dd));
Ydd_Phase2 = angle(Z13dd)/2/pi*360; 
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500    0.3250    0.0980];
figure(5)
clf
subplot(2,1,1);
%semilogx(f_s,Ydd_mag_db,'linewidth',1,'color',color_blue); grid on; hold on;
semilogx(f_s,Ydd_mag_db,'rx','linewidth',2,'color',color1,'MarkerSize',7); grid on; hold on; %grid off;
semilogx(f_s,Ydd_mag_db2,'rx','linewidth',2,'color',color2,'MarkerSize',7); grid on; hold on; %grid off;
ylabel('Magnitude (dB)')
title('$Y^{sys}_{dd}$','interpreter','latex','FontSize',12)
subplot(2,1,2);
%semilogx(f_s,Ydd_Phase,'linewidth',1,'color',color_blue); grid on; hold on;
semilogx(f_s,Ydd_Phase,'rx','linewidth',2,'color',color1,'MarkerSize',7); grid on; hold on; %grid off;
semilogx(f_s,Ydd_Phase2,'rx','linewidth',2,'color',color2,'MarkerSize',7); grid on; hold on; %grid off;
ylabel('Phase (degree)')
xlabel('Frequency (Hz)')