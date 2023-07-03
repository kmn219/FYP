% Calculate the admittance sensitivity matrix and nodal admittance matrix at the k-th labmda, 
% Final results will be numerical matrices.
% SensMat: eigenvalue-sensitivity matrix
% Ybus_val: network admittance value(lines + passive loads)
% Ynodal_val: value of the system nodal admittance matrix
% Yre_val: value of a rearranged nodal admittance matrix: diag->node, off
%          diag -> branch.
% Author: Yue Zhu

function [SensMat,Ybus_val,Ynodal_val,Yre_val] = SensitivityCal(ZminSS,Ek,lambda_rad,YbusObj)

%GmObj_Cell=evalin('base', 'GmObj_Cell');
%YbusObj=evalin('base', 'YbusObj');
N_Apparatus=evalin('base', 'N_Apparatus');
N_Bus=evalin('base', 'N_Bus');
GmDSS_Cell=evalin('base', 'GmDSS_Cell');
ApparatusBus=evalin('base', 'ApparatusBus');

%ZminSS = SimplusGT.WholeSysZ_cal(GmObj_Cell,YbusObj,N_Apparatus,N_Bus);

A=ZminSS.A;
B=ZminSS.B;
C=ZminSS.C;
[Phi,D]=eig(A);
Psi=inv(Phi); 
Mode=diag(D);

%% Admittance Sensitivity Matrix = -1* residue matrix
%SensMat_exp = zeros(N_Bus*2);
% for i=1:N_Bus*2
%     for j=1:N_Bus*2
%         SensMat_exp(i,j)=-1*C(i,:)*Phi(:,Ek)*Psi(Ek,:)*B(:,j);
%     end
% end
SensMat_exp = -1*C*Phi(:,Ek)*Psi(Ek,:)*B; % Residue matrix in expansion form

% pack up into d-q format
for i = 1:N_Bus
    for j = 1:N_Bus
        SensMat(i,j).dd = SensMat_exp(2*i-1,2*j-1);
        SensMat(i,j).dq = SensMat_exp(2*i-1,2*j);
        SensMat(i,j).qd = SensMat_exp(2*i,2*j-1);
        SensMat(i,j).qq = SensMat_exp(2*i,2*j);
    end
end

%% Node and branch Admittance values
[~,YbusDSS]=YbusObj.GetDSS(YbusObj);
Ybus_val_exp = evalfr(YbusDSS, lambda_rad); % Nodal admittance matrix Value (only passive part)

Ynodal_val_exp = Ybus_val_exp; % Nodal admittance matrix (include apparatus admittance)

for i=1:N_Apparatus
    bus_i = ApparatusBus{i};
    Ynodal_val_exp(2*bus_i-1,2*bus_i-1) = Ybus_val_exp(2*bus_i-1,2*bus_i-1)+ evalfr(GmDSS_Cell{i}(1,1),lambda_rad); %dd
    Ynodal_val_exp(2*bus_i-1,2*bus_i  ) = Ybus_val_exp(2*bus_i-1,2*bus_i)  + evalfr(GmDSS_Cell{i}(1,2),lambda_rad); %dq
    Ynodal_val_exp(2*bus_i  ,2*bus_i-1) = Ybus_val_exp(2*bus_i,  2*bus_i-1)+ evalfr(GmDSS_Cell{i}(2,1),lambda_rad); %qd
    Ynodal_val_exp(2*bus_i  ,2*bus_i  ) = Ybus_val_exp(2*bus_i,  2*bus_i)  + evalfr(GmDSS_Cell{i}(2,2),lambda_rad); %qq    
end

for i = 1:N_Bus
    for j = 1:N_Bus
        Ybus_val(i,j).dd = Ybus_val_exp(2*i-1,2*j-1);
        Ybus_val(i,j).dq = Ybus_val_exp(2*i-1,2*j);
        Ybus_val(i,j).qd = Ybus_val_exp(2*i,2*j-1);
        Ybus_val(i,j).qq = Ybus_val_exp(2*i,2*j);
        
        Ynodal_val(i,j).dd = Ynodal_val_exp(2*i-1,2*j-1);
        Ynodal_val(i,j).dq = Ynodal_val_exp(2*i-1,2*j);
        Ynodal_val(i,j).qd = Ynodal_val_exp(2*i,2*j-1);
        Ynodal_val(i,j).qq = Ynodal_val_exp(2*i,2*j);        
    end
end

% Rearranged network admittance matrix
for i = 1:N_Bus
    for j=1:N_Bus
        if i == j % node element : sum of a row / column
            Yre_val(i,i).dd = 0;%Ynodal_val(i,i).dd;
            Yre_val(i,i).dq = 0;%Ynodal_val(i,i).dq;
            Yre_val(i,i).qd = 0;%Ynodal_val(i,i).qd;
            Yre_val(i,i).qq = 0;%Ynodal_val(i,i).qq;
            for k = 1:N_Bus
               Yre_val(i,i).dd = Yre_val(i,i).dd+Ynodal_val(i,k).dd;
               Yre_val(i,i).dq = Yre_val(i,i).dq+Ynodal_val(i,k).dq;
               Yre_val(i,i).qd = Yre_val(i,i).qd+Ynodal_val(i,k).qd;
               Yre_val(i,i).qq = Yre_val(i,i).qq+Ynodal_val(i,k).qq;
            end
        else % branch element : inverse of the Ynodal
            Yre_val(i,j).dd = -1*Ynodal_val(i,j).dd;
            Yre_val(i,j).dq = -1*Ynodal_val(i,j).dq;
            Yre_val(i,j).qd = -1*Ynodal_val(i,j).qd;
            Yre_val(i,j).qq = -1*Ynodal_val(i,j).qq;
        end
    end
end



end