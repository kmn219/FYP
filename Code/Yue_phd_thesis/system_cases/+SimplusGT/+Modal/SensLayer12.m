function [SensLayer1, SensLayer2, Layer12] = SensLayer12(SensMatrix,Yre_val,modei,ModeSel)
N_Bus=evalin('base', 'N_Bus');
for i=1:N_Bus
    for j=1:N_Bus
        if i == j % node    
            SensLayer1(i,j)= SimplusGT.Frobenius_norm_dq(SensMatrix(i,j)) ...
                * SimplusGT.Frobenius_norm_dq(Yre_val(i,j));
            SensLayer2(i,j) = SimplusGT.inner_product_dq(SensMatrix(i,j),Yre_val(i,j));
        else % branch
            Sens.dd = SensMatrix(i,i).dd + SensMatrix(j,j).dd - SensMatrix(i,j).dd - SensMatrix(j,i).dd;
            Sens.dq = SensMatrix(i,i).dq + SensMatrix(j,j).dq - SensMatrix(i,j).dq - SensMatrix(j,i).dq;
            Sens.qd = SensMatrix(i,i).qd + SensMatrix(j,j).qd - SensMatrix(i,j).qd - SensMatrix(j,i).qd;
            Sens.qq = SensMatrix(i,i).qq + SensMatrix(j,j).qq - SensMatrix(i,j).qq - SensMatrix(j,i).qq;
            SensLayer1(i,j)= SimplusGT.Frobenius_norm_dq(Sens) * SimplusGT.Frobenius_norm_dq(Yre_val(i,j));
            SensLayer2(i,j) = SimplusGT.inner_product_dq(Sens,Yre_val(i,j));
        end
    end
end


%%
%%diagrams drawing
close(figure(2110+modei));
h=figure(2110+modei);
set(h,'position',[254.6,151.4,976.8,556]);

% Count=0;
% Layer2 results abs sum, for normalization.
Layer1Sum=0;
Layer2RealSum=0;
Layer2ImagSum=0;

nodeL1val = diag(SensLayer1);
Layer1Sum = (sum(SensLayer1,'all')+sum(nodeL1val,'all')) / 2;
nodeL2val = diag(SensLayer2);
Layer2RealSum = ( sum(abs(real(SensLayer2)),'all')+sum(abs(real(nodeL2val)),'all') )/2;
Layer2ImagSum = ( sum(abs(imag(SensLayer2)),'all')+sum(abs(imag(nodeL2val)),'all') )/2;


Count=0;
for k = 1:N_Bus % node    
       Count = Count + 1;
       Layer1(Count) = SensLayer1(k,k);
       Layer2.real(Count) = real(SensLayer2(k,k));
       Layer2.imag(Count) = imag(SensLayer2(k,k));
       Layer2.real_pu(Count) = real(SensLayer2(k,k))/Layer2RealSum;
       Layer2.imag_pu(Count) = imag(SensLayer2(k,k))/Layer2ImagSum;
       VecLegend{Count} = ['Node',num2str(k)];
       c(Count) = categorical({['Node',num2str(k)]});
       Layer12(Count).Component=['Node',num2str(k)];
       Layer12(Count).L1val_norm=Layer1(Count)/Layer1Sum;
       Layer12(Count).L2val_real_norm=Layer2.real_pu(Count);
       Layer12(Count).L2val_imag_norm=Layer2.imag_pu(Count);
end
for k = 1:N_Bus-1 % branch
    for j = k:N_Bus
        if k~=j && SensLayer1(k,j)~=0
           Count = Count + 1;
           Layer1(Count) = SensLayer1(k,j);
           Layer2.real(Count) = real(SensLayer2(k,j));
           Layer2.imag(Count) = imag(SensLayer2(k,j));
           Layer2.real_pu(Count) = real(SensLayer2(k,j))/Layer2RealSum;
           Layer2.imag_pu(Count) = imag(SensLayer2(k,j))/Layer2ImagSum;
           VecLegend{Count} = ['Branch',num2str(k),'-',num2str(j)];
           c(Count) = categorical({['Branch',num2str(k),'-',num2str(j)]});
           Layer12(Count).Component=['Branch',num2str(k),'-',num2str(j)];
           Layer12(Count).L1val_norm=Layer1(Count)/Layer1Sum;
           Layer12(Count).L2val_real_norm=Layer2.real_pu(Count);
           Layer12(Count).L2val_imag_norm=Layer2.imag_pu(Count);
        end
    end
end

Layer1 = Layer1/sum(Layer1);

clear title
subplot(2,2,[1,3]);
pie(Layer1,VecLegend);
title ('Admittance Sensitivity Level-1');
legend(VecLegend,'Location',[0.0486,0.121,0.1038,0.2437]);
colormap colorcube

subplot(2,2,2);
b=bar(c, Layer2.real_pu);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',-1:0.5:1);
%set(gca,'XTickLabel',[]);
title ('Admittance Sensitivity Level-2 Real (Normalized to 1)');
grid on;
% for i=1:Count
%     text(i-0.4,Layer2.real(i),num2str(Layer2.real(i)));
% end

subplot(2,2,4);
b=bar(c, Layer2.imag_pu);
set(gca,'YLim',[-1,1])
set(gca,'YTick',-1:0.5:1);
b.FaceColor = 'flat';
b.CData = [1,0.5,0];
title ('Admittance Sensitivity Level-2 Imaginary (Normalized to 1)');
grid on;
% for i=1:Count
%     text(i-0.4,Layer2.imag(i),num2str(Layer2.imag(i)));
% end

TitStr = ['Mode: ',num2str(ModeSel,'%.2f'), ' Hz'];
SimplusGT.mtit(TitStr, 'color',[1 0 0], 'xoff', -0.3);

end
