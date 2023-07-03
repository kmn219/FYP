% Calculate parameter paraticipatin factor of the selected devices, and all
% branches and passive loads.

function [Layer3_app,Layer3_bus] = SensLayer3(SensMatrix,Mode_rad,ApparatusSelL3All,Line_sel)

GmDSS_Cell=evalin('base', 'GmDSS_Cell');
Para = evalin('base', 'Para');
Ts = evalin('base', 'Ts');
ListBus=evalin('base', 'ListBus');
ApparatusPowerFlow= evalin('base', 'ApparatusPowerFlow');
ApparatusType=evalin('base', 'ApparatusType');
ApparatusBus=evalin('base', 'ApparatusBus');

%if length(ApparatusSelL3All)==0
%    error('you need to choose an Apparatus in Sheet-2 colum L in Modalconfig excel file')
    
for AppCount = 1:length(ApparatusSelL3All)
    AppSel = ApparatusSelL3All(AppCount);
    1i;
    YmValOrig.dd=evalfr(GmDSS_Cell{AppSel}(1,1),Mode_rad);
    YmValOrig.dq=evalfr(GmDSS_Cell{AppSel}(1,2),Mode_rad);
    YmValOrig.qd=evalfr(GmDSS_Cell{AppSel}(2,1),Mode_rad);
    YmValOrig.qq=evalfr(GmDSS_Cell{AppSel}(2,2),Mode_rad);
    
    ParamName = fieldnames(Para{AppSel});
    %Residue_ = Residue(AppSel);
    %perturb the parameters one by one.
    for k=1:length(ParamName)
        ParaNew = Para;
        ParaSel = getfield(Para{AppSel},ParamName{k}); % extract the parameter
        delta_para = 1e-5*(1+abs(ParaSel));
        ParaPerturb = ParaSel + delta_para ; % add perturabation
        ParaNew = setfield(ParaNew{AppSel}, ParamName{k}, ParaPerturb); % update the parameter  
   
        [~,GmDSS_Cell_New,~,~,~,~,~,~,~] ...    % get the new parameter
        = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{AppSel},ApparatusType{AppSel},...
                            ApparatusPowerFlow{AppSel},ParaNew,Ts,ListBus);                        
        1i;
        YmValNew.dd=evalfr(GmDSS_Cell_New(1,1),Mode_rad);
        YmValNew.dq=evalfr(GmDSS_Cell_New(1,2),Mode_rad);
        YmValNew.qd=evalfr(GmDSS_Cell_New(2,1),Mode_rad);
        YmValNew.qq=evalfr(GmDSS_Cell_New(2,2),Mode_rad);                                                
        
        Layer3_app(AppCount).Apparatus={['Apparatus',num2str(AppSel)]};
        Layer3_app(AppCount).Result(k).ParaName = ParamName(k);
        Layer3_app(AppCount).Result(k).DeltaY.dd = (YmValNew.dd - YmValOrig.dd)/(delta_para);
        Layer3_app(AppCount).Result(k).DeltaY.dq = (YmValNew.dq - YmValOrig.dq)/(delta_para);
        Layer3_app(AppCount).Result(k).DeltaY.qd = (YmValNew.qd - YmValOrig.qd)/(delta_para);
        Layer3_app(AppCount).Result(k).DeltaY.qq = (YmValNew.qq - YmValOrig.qq)/(delta_para);
        DeltaY = Layer3_app(AppCount).Result(k).DeltaY;
        
        Layer3_app(AppCount).Result(k).DLambda_rad = SimplusGT.inner_product_dq(DeltaY, SensMatrix(AppSel,AppSel));
        DLambda_Hz=Layer3_app(AppCount).Result(k).DLambda_rad/(2*pi);
        Layer3_app(AppCount).Result(k).DLambdaRho_Hz=DLambda_Hz;
        Layer3_app(AppCount).Result(k).DLambdaRho_pu_Hz=DLambda_Hz*ParaSel;
    end
end


%% load and branch
% ZminSS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i, Port_v);
% [SensMatrix, Ybus_val, Ynodal_val, Yre_val] = SimplusGT.Modal.SensitivityCal(ZminSS,Ek,mode_rad);

ListLineNew=evalin('base', 'ListLineNew');
ListBus=evalin('base', 'ListBus');
Wbase=evalin('base', 'Wbase');
YbusDSS=evalin('base', 'YbusDSS');
N_Bus=evalin('base', 'N_Bus');

FB    = ListLineNew(:,1);              % From bus number...
TB    = ListLineNew(:,2);              % To bus number...
% Rlist = ListLineNew(:,3);              % Resistance,  R...
% Xlist = ListLineNew(:,4);              % Inductance,  wL...
% Blist = ListLineNew(:,5);              % Capacitance, wC...
% Glist = ListLineNew(:,6);              % Conductance, G...
% Tlist = ListLineNew(:,7);              % Turns ratio, T

count=1;
[List_row,~] = size(ListLineNew);
CompoName=['R','X','B','G'];

for m = 1:length(Line_sel) % all rows
    k = Line_sel(m);
    for j =1: 4 % R,X,B,G
        if ListLineNew(k,j+2)~=inf && ListLineNew(k,j+2)~=0 % not inf, not 0
            ListLineNew_New = ListLineNew;
            Zorig = ListLineNew(k,j+2);
            DZ = 1e-5*(1+abs(Zorig));
            Znew=Zorig+DZ;
            ListLineNew_New(k,j+2)=ListLineNew(k,j+2)+DZ;
            if FB(k)==TB(k)
                name_z = ['shunt',num2str(FB(k)),'.',CompoName(j)];
            else
                name_z = ['branch',num2str(FB(k)),'-',num2str(TB(k)),'.',CompoName(j)];
            end
            
            %%%%%%%%%%
            Sens_branch.dd=SensMatrix(FB(k),FB(k)).dd+SensMatrix(TB(k),TB(k)).dd-...
                SensMatrix(FB(k),TB(k)).dd-SensMatrix(TB(k),FB(k)).dd;
            
            Sens_branch.dq=SensMatrix(FB(k),FB(k)).dq+SensMatrix(TB(k),TB(k)).dq-...
                SensMatrix(FB(k),TB(k)).dq-SensMatrix(TB(k),FB(k)).dq;
            
            Sens_branch.qd=SensMatrix(FB(k),FB(k)).qd+SensMatrix(TB(k),TB(k)).qd-...
                SensMatrix(FB(k),TB(k)).qd-SensMatrix(TB(k),FB(k)).qd;
            
            Sens_branch.qq=SensMatrix(FB(k),FB(k)).qq+SensMatrix(TB(k),TB(k)).qq-...
                SensMatrix(FB(k),TB(k)).qq-SensMatrix(TB(k),FB(k)).qq;
            if j==1
                Dzx=[Znew-Zorig,0; 0, Znew-Zorig]/DZ;
                Dy_=inv(Dzx);
                Dy.dd=Dy_(1,1);
                Dy.dq=Dy_(1,2);
                Dy.qd=Dy_(2,1);
                Dy.qq=Dy_(2,2);
            elseif j==2
                Dzx=[0, 1i*(Znew-Zorig); -1i*(Znew-Zorig), 0]/DZ;
                Dy_=inv(Dzx);
                Dy.dd=Dy_(1,1);
                Dy.dq=Dy_(1,2);
                Dy.qd=Dy_(2,1);
                Dy.qq=Dy_(2,2);
            else
                error('this line is not support for Layer-3 for now, choose another line in ModalAnalysisExe.m!')
            end
            D_lambda_rad = SimplusGT.inner_product_dq(Sens_branch,Dy);
            %%%%%%%%%%%%%%%%%
            
%             Delta_Y_exp=Delta_Y_exp/DZ;
%             [~,YbusDSS_new,~] = SimplusGT.Toolbox.YbusCalcDss(ListBus, ListLineNew_New, Wbase);
%             YbusNew_val_exp = evalfr(YbusDSS_new, Mode_rad);
%             YbusOrig_val_exp = evalfr(YbusDSS, Mode_rad);
% %             ZminSS_New = SimplusGT.WholeSysZ_cal(GmObj,YbusObjNew,Port_i, Port_v);
% %             [~, ~, Ynodal_val_new,~]=SimplusGT.Modal.SensitivityCal(ZminSS_New,1,Mode_rad,YbusObjNew);
%             Delta_Y_exp = (YbusNew_val_exp-YbusOrig_val_exp)/DZ;
%             for i = 1:length(Delta_Y_exp)/2
%                 for j = 1:length(Delta_Y_exp)/2
%                     Delta_Y(i,j).dd = Delta_Y_exp(2*i-1,2*j-1);
%                     Delta_Y(i,j).dq = Delta_Y_exp(2*i-1,2*j);
%                     Delta_Y(i,j).qd = Delta_Y_exp(2*i,2*j-1);
%                     Delta_Y(i,j).qq = Delta_Y_exp(2*i,2*j); 
%                 end
%             end
        
            %SensMatrix
            %Delta_Y = SimplusGT.dqStrutMinus(YbusNew_val,YbusOrig_val);
%            D_lambda_rad = SimplusGT.TraceFo_dq(SensMatrix,Delta_Y); % trace calculation.
                        
            D_lambda_Hz_pu = D_lambda_rad * 2*pi * Zorig ;
            Layer3_bus(count).component = name_z;
            Layer3_bus(count).D_lambda_rad = D_lambda_rad;
            Layer3_bus(count).D_lambda_Hz_pu = D_lambda_Hz_pu;
            count = count+1;
        end
    end
end

%[YbusObj,YbusDSS,~] = SimplusGT.Toolbox.YbusCalcDss(ListBus,ListLineNew,Wbase);



% D_Ynodal = Ynodal_val_new - Ynodal_val_orig;
% D_lambda = SimplusGT.TraceFo_dq(SensMatrix,D_Ynodal);

end
