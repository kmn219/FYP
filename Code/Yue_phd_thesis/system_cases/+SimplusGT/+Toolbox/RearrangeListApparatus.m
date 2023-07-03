% This function re-arranges the netlist data of apparatuses.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The apparatus model is in load convention.
%
% The apparautus bus index, type index, and parameters are saved in three
% different cells. The reason of doing this is that, the number of
% apparatuses might not be same to the number of buses, when the power
% system has interlink apparatuses. In this case, the interlink apparatus
% will have two bus indice.
%
% "app" means apparatus.

function [AppBusCell,AppTypeCell,ParaCell,N_App] = RearrangeListApparatus(UserData,W0,ListBus)

%% New excel in xlsm format
% Check if the xlsm format is used
if strcmp(UserData,'UserData.xlsm')
    NewExcel = 1;
else
    NewExcel = 0;
end

%% Load data
[ListApp,ListAppChar]	 = xlsread(UserData,'Apparatus');

%% Prepare
if NewExcel == 1
    ListApp = ListApp(3:end,:);       % Remove the first two lines
    ListApp = ListApp(:,[1:2,4:end]); % Remove the 3rd colomn of subtype
    [rmax,cmax] = size(ListApp);
    ListAppNew = NaN(rmax/2,cmax);
    for r = 1:(rmax/2)
        ListAppNew(r,1) = ListApp(2*r-1,1);         % Move bus number
        ListAppNew(r,2) = ListApp(2*r-1,2);         % Move bus number
        ListAppNew(r,3:end) = ListApp(2*r,3:end);   % Move others
    end
   ListApp = ListAppNew;            % Update
end
ListAppBusChar = ListAppChar(:,1);
for k1 = 1:length(ListAppBusChar)
    if strcmpi(ListAppBusChar{k1},'Bus No.')
        break;
    end
end
ListAppBusChar = ListAppBusChar(k1+1:end);
if NewExcel == 1
    ListAppBusChar = ListAppBusChar(1:2:end,:);
end

%% Rearrange data
% Notes:
% This section will convert the ListAppBus and ListAppType to AppBusCell
% and AppTypeCell, i.e., array type to cell type.

[N_App,ColumnMax_Apparatus] = size(ListApp);
ListAppBus = ListApp(:,1);
ListAppType = ListApp(:,2);

% Get the apparatus bus in cell form
for n = 1:N_App
    if ~isnan(ListAppBus(n))  % If not NaN, then the bus is a scalar rather han an array in char type
        AppBusCell{n} = ListAppBus(n);
    else
        AppBusCell{n} = str2num(ListAppBusChar{n});
        [~,~,AreaType]= SimplusGT.Toolbox.CheckBus(AppBusCell{n}(1),ListBus);
        
        % Notes:
        % If the first bus is dc bus, then swap. This ensures the first one
        % is ac bus.
        if AreaType == 2 
            [AppBusCell{n}(1),AppBusCell{n}(2)] = deal(AppBusCell{n}(2),AppBusCell{n}(1));
        end
    end
end

% Get the apparatus type in cell form
for n = 1:N_App
    AppTypeCell{n} = ListApp(n,2);
end

% Add floating bus
N_Bus = length(ListBus(:,1));
for m = 1:N_Bus
    BusIndex = ListBus(m,1);
    ExistApparatus = SimplusGT.CellFind(AppBusCell,BusIndex);
    if isempty(ExistApparatus)
        % The bus has no apparatus, i.e., an ac or dc floating bus.
        N_App = N_App+1;
        AppBusCell{N_App} = BusIndex;   % Add a new bus index
        [~,~,AreaType] = SimplusGT.Toolbox.CheckBus(BusIndex,ListBus);  % Check if an ac or dc bus
        if AreaType == 1
            AppTypeCell{N_App} = 100;       % Ac floating bus
        elseif AreaType == 2
            AppTypeCell{N_App} = 1100;      % Dc floating bus
        else
            error(['Error: Error AreaType.']);
        end
    else
        % The bus has an apparatus already, no need to do anything then.
    end
end

% Error check
if (ColumnMax_Apparatus>50)
    error(['Error: Apparatus data overflow.']); 
end

[~,ModeBus] = SimplusGT.CellMode(AppBusCell);
if ModeBus~=1
    error(['Error: For each bus, one and only one apparatus has to be connected.']);
end

%% Default AC apparatus data
% ======================================
% Synchronous generator
% ======================================
Para0000.J  = 3.5;
Para0000.D  = 1;
Para0000.wL = 0.05;
Para0000.R  = 0.01;
Para0000.w0 = W0;

% ======================================
% Grid-following VSI (PLL-controlled)
% ======================================
% Dc link
Para0010.V_dc       = 2.5;
Para0010.C_dc       = 1.25;         % 2*0.1*Para0010.V_dc^2;
Para0010.f_v_dc     = 5;            % (Hz) bandwidth, vdc

% Ac filter
Para0010.wLf        = 0.03;
Para0010.R          = 0.01;

% PLL
Para0010.f_pll      = 5;            % (Hz) bandwidth, PLL
Para0010.f_tau_pll  = 300;          % (Hz) bandwidth, PLL low pass filter

% Current loop
Para0010.f_i_dq     = 600;      	% (Hz) bandwidth, idq
Para0010.w0         = W0;   

% ======================================
% Grid-forming VSI (Droop-Controlled)
% ======================================
Para0020.wLf    =0.05;
Para0020.Rf     =0.05/5;
Para0020.wCf    =0.02;
Para0020.wLc    =0.01;
Para0020.Rc     =0.01/5;
Para0020.Xov    =0.01;
Para0020.Dw     =0.05;
Para0020.fdroop =5;    % (Hz) droop control bandwidth
Para0020.fvdq   =300;   % (Hz) vdc bandwidth
Para0020.fidq   =600;   % current control bandwidth
Para0020.w0     = W0;

% ======================================
% Synchronous generato ----- Full Model
% ======================================
Para0030.X=0.0125;
Para0030.R=0;
Para0030.Xd=0.1; %synchronous reactance in d axis
Para0030.Xd1=0.031; %transient reactance
Para0030.Xd2=0.025; %subtransient reactance
Para0030.Td1=10.2; %d-axis open circuit transient time constant
Para0030.Td2=0.05; %d-axis open circuit sub-transient time constant
Para0030.Xq=0.069;
Para0030.Xq1=0.0416667;
Para0030.Xq2=0.025;
Para0030.Tq1=1.5;
Para0030.Tq2=0.035;
Para0030.H=42;
Para0030.D=0;
Para0030.TR=0.01;
Para0030.KA=1;
Para0030.TA=0.02;
Para0030.VRmax=10;
Para0030.VRmin=-10;
Para0030.KE=1;
Para0030.TE=0.785;
Para0030.E1=3.9267;
Para0030.SE1=0.07;
Para0030.E2=5.2356;
Para0030.SE2=0.91;
Para0030.KF=0.03;
Para0030.TF=1;
Para0030.KP=200;
Para0030.KI=50;
Para0030.KD=50;
Para0030.TD=0.01;
Para0030.KPSS=20;
Para0030.TW=15;
Para0030.T11=0.15;
Para0030.T12=0.04;
Para0030.T21=0.15;
Para0030.T22=0.04;
Para0030.T31=0.15;
Para0030.T32=0.04;
Para0030.VSSmax=0.2;
Para0030.VSSmin=-0.05;
Para0030.Rgov=0.05;
Para0030.T1gov=0.8;
Para0030.T2gov=1.8;
Para0030.T3gov=6;
Para0030.Dtgov=0.2;

% ======================================
% Synchronous generato ----- Full Model ---- Cyprus Model
% ======================================
Para0040.Sbase_SM=1;
Para0040.X=0.0125;
Para0040.R=0;
Para0040.Xd=0.1; %synchronous reactance in d axis
Para0040.Xd1=0.031; %transient reactance
Para0040.Xd2=0.025; %subtransient reactance
Para0040.Td1=10.2; %d-axis open circuit transient time constant
Para0040.Td2=0.05; %d-axis open circuit sub-transient time constant
Para0040.Xq=0.069;
Para0040.Xq1=0.0416667;
Para0040.Xq2=0.025;
Para0040.Tq1=1.5;
Para0040.Tq2=0.035;
Para0040.H=42;
Para0040.D=0;
Para0040.Dpu=0;
Para0040.SG10=1;
Para0040.SG12=1;
%AVR
Para0040.Tr=0;
Para0040.Ka=50;
Para0040.Ta=0.06;
Para0040.Vrmax=1;
Para0040.Vrmin=-1;
Para0040.Ke=-0.0465;
Para0040.Te=0.520;
Para0040.Kf=0.0832;
Para0040.Tf=1.0;
Para0040.E1=3.24;
Para0040.SEE1=0.072;
Para0040.E2=4.320;
Para0040.SEE2=0.2821;
%PSS
Para0040.T1=5;
Para0040.T2=0.6;
Para0040.T3=3;
Para0040.T4=0.5;
Para0040.Tw=10;
Para0040.Kpss=1;
Para0040.Vpssmin=-0.2;
Para0040.Vpssmax=0.2;
%Governor
Para0040.Rgov=0.011;
Para0040.T1gov=0.1;
Para0040.T2gov=0;
Para0040.T3gov=0.3;
Para0040.T4gov=0.05;
Para0040.T5gov=10;
Para0040.Fgov=0.25;
Para0040.Pmax_gov=1;
Para0040.w0 = W0;


% ======================================
% Ac infinite bus (short-circuit in small-signal)
% ======================================
Para0090 = [];

% ======================================
% Ac floating bus (open-circuit)
% ======================================
Para0100 = [];

%% Default DC apparatus data
% ======================================
% Grid-feeding buck
% ======================================
Para1010.Vdc  = 2;
Para1010.Cdc  = 0.8;
Para1010.wL   = 0.05;
Para1010.R    = 0.05/5;
Para1010.fi   = 600;
Para1010.fvdc = 5;
Para1010.w0   = W0;

% ======================================
% Dc infinite bus (short-circuit in small-signal)
% ======================================
Para1090 = [];

% ======================================
% Dc floating bus (open-circuit)
% ======================================
Para1100 = [];

%% Default hybrid apparatus data
% ======================================
% Interlink ac-dc converter
% ======================================
Para2000.C_dc   = 1.6;
Para2000.wL_ac  = 0.05;
Para2000.R_ac   = 0.01;
Para2000.wL_dc  = 0.02;
Para2000.R_dc   = 0.02/5;
Para2000.fidq   = 600;
Para2000.fvdc   = 5;
Para2000.fpll   = 5;
Para2000.w0     = W0;   

%% Re-arrange apparatus data

% Find the index of user-defined data
ListApparatus_NaN = isnan(ListApp(:,3:ColumnMax_Apparatus));    % Find NaN
[row,column] = find(ListApparatus_NaN == 0);     
column = column+2;

% Initialize the apparatus parameters by default parameters
for i = 1:N_App
    AppBus   = AppBusCell{i};
    AppType  = AppTypeCell{i};
    switch floor(AppType/10)
        % ### AC apparatuses
        case 0     
            ParaCell{i} = Para0000;     % Synchronous machine
        case 1
            ParaCell{i} = Para0010;     % Grid-following inverter
      	case 2
            ParaCell{i} = Para0020;     % Grid-forming inverter
        case 3
            ParaCell{i} = Para0030;     % Synchronous Machine full model (68-bus model)
        case 4
            ParaCell{i} = Para0040;     % Synchronous Machine full model (Cyprus model)
            % Yue's Full-Order Machine
        case 9
            ParaCell{i} = Para0090;     % Ac inifnite bus
        case 10
            ParaCell{i} = Para0100;     % Ac floating bus, i.e., no apparatus
        
        % ### DC apparatuses
        case 101
            ParaCell{i} = Para1010;     % Grid-following buck
        case 109
            ParaCell{i} = Para1090;     % Dc infinite bus
        case 110
            ParaCell{i} = Para1100;     % Ac floating bus, i.e., no apparatus
            
        % ### Hybrid ac-dc apparatuses
        case 200
            ParaCell{i} = Para2000;     % Interlinking ac-dc converter
            
        % ### Error check
        otherwise
            error(['Error: apparatus type, bus ' num2str(AppBus) ' type ' num2str(AppType) '.']);
    end
end

% Replace the default data by customized data
%
% Notes: 
% This method can reduce the calculation time of "for" loop. The "for" loop
% runs only when "row" is not empty.
%
% The sequence of cases are determined by the excel form. This also
% decouples the sequence between the excel form and the system object.
for i = 1:length(row)
  	AppBus   = AppBusCell{row(i)};
	AppType	= AppTypeCell{row(i)};
 	UserValue 	= ListApp(row(i),column(i));     % Customized value
    SwitchFlag = column(i)-2;                   	% Find the updated parameter
  	if floor(AppType/10) == 0                    % Synchronous machine
        switch SwitchFlag 
         	case 1; ParaCell{row(i)}.J  = UserValue;
            case 2; ParaCell{row(i)}.D  = UserValue;
            case 3; ParaCell{row(i)}.wL = UserValue;
            case 4; ParaCell{row(i)}.R  = UserValue; 
            otherwise
                error(['Error: paramter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif (floor(AppType/10) == 1)              % Grid-following inverter
        switch SwitchFlag
            case 1; ParaCell{row(i)}.V_dc   = UserValue;
            case 2; ParaCell{row(i)}.C_dc   = UserValue;
            case 3; ParaCell{row(i)}.wLf     = UserValue;
            case 4; ParaCell{row(i)}.R      = UserValue;
            case 5; ParaCell{row(i)}.f_v_dc = UserValue;
            case 6; ParaCell{row(i)}.f_pll  = UserValue;
            case 7; ParaCell{row(i)}.f_i_dq = UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 2                % Grid-forming inverter
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.wLf     = UserValue;
          	case 2;  ParaCell{row(i)}.Rf      = UserValue;
          	case 3;  ParaCell{row(i)}.wCf     = UserValue;
           	case 4;  ParaCell{row(i)}.wLc  	  = UserValue;
         	case 5;  ParaCell{row(i)}.Rc  	  = UserValue;
           	case 6;  ParaCell{row(i)}.Xov 	  = UserValue;
            case 7;  ParaCell{row(i)}.Dw      = UserValue;
            case 8;  ParaCell{row(i)}.fdroop  = UserValue;
          	case 9;  ParaCell{row(i)}.fvdq    = UserValue;
          	case 10; ParaCell{row(i)}.fidq    = UserValue; 
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
        
    elseif floor(AppType/10) == 3 %full model
        switch SwitchFlag
            case 1; ParaCell{row(i)}.X=UserValue;
            case 2; ParaCell{row(i)}.R=UserValue;
            case 3; ParaCell{row(i)}.Xd=UserValue; 
            case 4; ParaCell{row(i)}.Xd1=UserValue; 
            case 5; ParaCell{row(i)}.Xd2=UserValue; 
            case 6; ParaCell{row(i)}.Td1=UserValue; 
            case 7; ParaCell{row(i)}.Td2=UserValue;
            case 8; ParaCell{row(i)}.Xq=UserValue;
            case 9; ParaCell{row(i)}.Xq1=UserValue;
            case 10;ParaCell{row(i)}.Xq2=UserValue;
            case 11;ParaCell{row(i)}.Tq1=UserValue;
            case 12;ParaCell{row(i)}.Tq2=UserValue;
            case 13;ParaCell{row(i)}.H=UserValue;
            case 14;ParaCell{row(i)}.D=UserValue;
            case 15;ParaCell{row(i)}.TR=UserValue;
            case 16;ParaCell{row(i)}.KA=UserValue;
            case 17;ParaCell{row(i)}.TA=UserValue;
            case 18;ParaCell{row(i)}.VRmax=UserValue;
            case 19;ParaCell{row(i)}.VRmin=UserValue;
            case 20;ParaCell{row(i)}.KE=UserValue;
            case 21;ParaCell{row(i)}.TE=UserValue;
            case 22;ParaCell{row(i)}.E1=UserValue;
            case 23;ParaCell{row(i)}.SE1=UserValue;
            case 24;ParaCell{row(i)}.E2=UserValue;
            case 25;ParaCell{row(i)}.SE2=UserValue;
            case 26;ParaCell{row(i)}.KF=UserValue;
            case 27;ParaCell{row(i)}.TF=UserValue;
            case 28;ParaCell{row(i)}.KP=UserValue;
            case 29;ParaCell{row(i)}.KI=UserValue;
            case 30;ParaCell{row(i)}.KD=UserValue;
            case 31;ParaCell{row(i)}.TD=UserValue;
            case 32;ParaCell{row(i)}.KPSS=UserValue;
            case 33;ParaCell{row(i)}.TW=UserValue;
            case 34;ParaCell{row(i)}.T11=UserValue;
            case 35;ParaCell{row(i)}.T12=UserValue;
            case 36;ParaCell{row(i)}.T21=UserValue;
            case 37;ParaCell{row(i)}.T22=UserValue;
            case 38;ParaCell{row(i)}.T31=UserValue;
            case 39;ParaCell{row(i)}.T32=UserValue;
            case 40;ParaCell{row(i)}.VSSmax=UserValue;
            case 41;ParaCell{row(i)}.VSSmin=UserValue;
            case 42;ParaCell{row(i)}.Rgov=UserValue;
            case 43;ParaCell{row(i)}.T1gov=UserValue;
            case 44;ParaCell{row(i)}.T2gov=UserValue;
            case 45;ParaCell{row(i)}.T3gov=UserValue;
            case 46;ParaCell{row(i)}.Dtgov=UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    
    elseif floor(AppType/10) == 4 %full model Cyprus
        switch SwitchFlag
            case 1; ParaCell{row(i)}.Sbase_SM=UserValue;
            case 2; ParaCell{row(i)}.X=UserValue;
            case 3; ParaCell{row(i)}.R=UserValue;
            case 4; ParaCell{row(i)}.Xd=UserValue; %synchronous reactance in d axis
            case 5; ParaCell{row(i)}.Xd1=UserValue; %transient reactance
            case 6; ParaCell{row(i)}.Xd2=UserValue; %subtransient reactance
            case 7; ParaCell{row(i)}.Td1=UserValue; %d-axis open circuit transient time constant
            case 8; ParaCell{row(i)}.Td2=UserValue; %d-axis open circuit sub-transient time constant
            case 9; ParaCell{row(i)}.Xq=UserValue;
            case 10; ParaCell{row(i)}.Xq1=UserValue;
            case 11; ParaCell{row(i)}.Xq2=UserValue;
            case 12; ParaCell{row(i)}.Tq1=UserValue;
            case 13; ParaCell{row(i)}.Tq2=UserValue;
            case 14; ParaCell{row(i)}.H=UserValue;
            case 15; ParaCell{row(i)}.D=UserValue;
            case 16; ParaCell{row(i)}.Dpu=UserValue;
            case 17; ParaCell{row(i)}.SG10=UserValue;
            case 18; ParaCell{row(i)}.SG12=UserValue; 
            %AVR
            case 19; ParaCell{row(i)}.Tr=UserValue;
            case 20; ParaCell{row(i)}.Ka=UserValue;
            case 21; ParaCell{row(i)}.Ta=UserValue;
            case 22; ParaCell{row(i)}.Vrmax=UserValue;
            case 23; ParaCell{row(i)}.Vrmin=UserValue;
            case 24; ParaCell{row(i)}.Ke=UserValue;
            case 25; ParaCell{row(i)}.Te=UserValue;
            case 26; ParaCell{row(i)}.Kf=UserValue;
            case 27; ParaCell{row(i)}.Tf=UserValue;
            case 28; ParaCell{row(i)}.E1=UserValue;
            case 29; ParaCell{row(i)}.SEE1=UserValue;
            case 30; ParaCell{row(i)}.E2=UserValue;
            case 31; ParaCell{row(i)}.SEE2=UserValue;
            %PSS
            case 32; ParaCell{row(i)}.T1=UserValue;
            case 33; ParaCell{row(i)}.T2=UserValue;
            case 34; ParaCell{row(i)}.T3=UserValue;
            case 35; ParaCell{row(i)}.T4=UserValue;
            case 36; ParaCell{row(i)}.Tw=UserValue;
            case 37; ParaCell{row(i)}.Kpss=UserValue;
            case 38; ParaCell{row(i)}.Vpssmin=UserValue;
            case 39; ParaCell{row(i)}.Vpssmax=UserValue;
            %Governor
            case 40; ParaCell{row(i)}.Rgov=UserValue;
            case 41; ParaCell{row(i)}.T1gov=UserValue;
            case 42; ParaCell{row(i)}.T2gov=UserValue;
            case 43; ParaCell{row(i)}.T3gov=UserValue;
            case 44; ParaCell{row(i)}.T4gov=UserValue;
            case 45; ParaCell{row(i)}.T5gov=UserValue;
            case 46; ParaCell{row(i)}.Fgov=UserValue;
            case 47; ParaCell{row(i)}.Pmax_gov=UserValue;
            %Vbase    
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
        
    elseif floor(AppType/10) == 101 % Grid-feeding buck
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.Vdc   = UserValue;
          	case 2;  ParaCell{row(i)}.Cdc   = UserValue;
          	case 3;  ParaCell{row(i)}.wLf   = UserValue;
           	case 4;  ParaCell{row(i)}.R  	= UserValue;
         	case 5;  ParaCell{row(i)}.fi  	= UserValue;
           	case 6;  ParaCell{row(i)}.fvdc 	= UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 200 % Interlink ac-dc converter
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.C_dc  = UserValue;
            case 2;  ParaCell{row(i)}.wL_ac = UserValue;
            case 3;  ParaCell{row(i)}.R_ac  = UserValue;
            case 4;  ParaCell{row(i)}.wL_dc = UserValue;
            case 5;  ParaCell{row(i)}.R_dc  = UserValue;
            case 6;  ParaCell{row(i)}.fidq  = UserValue;
            case 7;  ParaCell{row(i)}.fvdc  = UserValue;
            case 8;  ParaCell{row(i)}.fpll  = UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    end
end

%% The re-order can only be done here
% For a single-area pure ac system, we re-order the apparatuses 
AreaTypeCheck = find(ListBus(:,12) == 2, 1);
AreaNoCheck = find(ListBus(:,11) == 2, 1);
if isempty(AreaTypeCheck) && isempty(AreaNoCheck)
    for j = 1:N_App
        [~,AppIndex] = SimplusGT.CellFind(AppBusCell,j);
        AppBusCellNew{j}    = AppBusCell{AppIndex};
        AppTypeCellNew{j}   = AppTypeCell{AppIndex};
        ParaCellNew{j}      = ParaCell{AppIndex};
    end
    AppBusCell  = AppBusCellNew;
    AppTypeCell = AppTypeCellNew;
    ParaCell    = ParaCellNew;
end
end