% This function adds apparatuses into simulink model.

% Author(s): Yitong Li

function [FullName_Apparatus,Name_Apparatus,Pos_Apparatus] = SimAddApparatus(Name_Model,Name_LibFile,Size_Apparatus,Shift_Apparatus,Pos_Bus,ApparatusBus,ApparatusType,ListAdvance)

% Organize data
DiscreMethod = ListAdvance(1);
LinearizationTimes = ListAdvance(2);
DiscreDampingFlag = ListAdvance(3);
DirectFeedthrough = ListAdvance(4);

N_Apparatus = length(ApparatusType);

% Add apparatus
for i = 1:N_Apparatus
    if ApparatusType{i}~=0100 && ApparatusType{i}~=1100
        
        % Get the bus index of apparatus
        Bus = ApparatusBus{i};
        
        switch floor(ApparatusType{i}/10)
            % ### Ac apparatus
            case 000
                Name_Apparatus{i} = ['SM' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Synchronous Machine (dq-Frame System Object)'],FullName_Apparatus{i});
            case 001
                Name_Apparatus{i} = ['VSI-PLL' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                if ApparatusType{i}<18
                    add_block([Name_LibFile '/Grid-Following Voltage-Source Inverter (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==18
                    add_block([Name_LibFile '/Grid-Following Voltage-Source Inverter (dq-Frame System Object)_old'],FullName_Apparatus{i});
                else
                    add_block([Name_LibFile '/Grid-Following Inverter (alpha/beta System Object)'],FullName_Apparatus{i});
                end
            case 002
                Name_Apparatus{i} = ['VSI-Droop' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                if ApparatusType{i}~=28
                    add_block([Name_LibFile '/Grid-Forming Voltage-Source Inverter (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==28
                    add_block([Name_LibFile '/Grid-Forming Voltage-Source Inverter (dq-Frame System Object)_old'],FullName_Apparatus{i});
                end
                
            case 003
                %Name_Apparatus{i} = ['SM_Full' num2str(Bus)];
                %FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                if ApparatusType{i}==30
                    Name_Apparatus{i} = ['SM_8th' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SM (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==31
                    Name_Apparatus{i} = ['SM_8th_AVR' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SMAVR (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==32
                    Name_Apparatus{i} = ['SM_8th_AVRPSS' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SMAVRPSS (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==33
                    Name_Apparatus{i} = ['SM_8th_GOV' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SM_GOV (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==34
                    Name_Apparatus{i} = ['SM_8th_AVRGOV' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SMAVRGOV (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==35
                    Name_Apparatus{i} = ['SM_8th_AVRPSSGOV' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SMAVRPSSGOV (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==38
                    Name_Apparatus{i} = ['SM_8th_FS' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SMGOV_FS (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==39
                    Name_Apparatus{i} = ['SM_8th_AVRPSSGOV_FS' num2str(i)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/SynchronousMachineFull_SMAVRPSSGOV_FS (dq-Frame System Object)'],FullName_Apparatus{i});
                end
                
            case 004 % Cyprus models
                Name_Apparatus{i} = ['SM' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                if ApparatusType{i}==41
                    add_block([Name_LibFile '/SynchronousMachineFull_Cyp_SMAVR (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==42
                    add_block([Name_LibFile '/SynchronousMachineFull_Cyp_SMAVRPSSGOV (dq-Frame System Object)'],FullName_Apparatus{i});
                end
            case 009
            	Name_Apparatus{i} = ['Inf-Bus' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/AC Infinite Bus (Voltage Type)'],FullName_Apparatus{i});
                
            % ### Dc apparatus
            case 101
                Name_Apparatus{i} = ['Buck' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Grid-Feeding Buck (System Object)'],FullName_Apparatus{i});
        	case 109
            	Name_Apparatus{i} = ['Inf-Bus' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/DC Infinite Bus (Voltage Type)'],FullName_Apparatus{i});
                
            % ### Interlink
            case 200
              	Name_Apparatus{i} = ['Interlink' num2str(Bus(1)) '-' num2str(Bus(2))];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Interlink AC-DC (System Object)'],FullName_Apparatus{i});
                
          	% ### Error check
            otherwise
                error(['Error: ApparatusType ' num2str(ApparatusType{i}) '.']);
        end
        
        % Set position
       	% The position of apparatus is set by referring to the position of correpsonding bus
        if 0<=ApparatusType{i} && ApparatusType{i}<=90    
            % For ac apparatuses
            Pos_Apparatus{i} = Pos_Bus{Bus} + Shift_Apparatus;
            set_param(FullName_Apparatus{i},'position',[Pos_Apparatus{i},Pos_Apparatus{i}+Size_Apparatus]);
        elseif 1000<=ApparatusType{i} && ApparatusType{i}<=1090   
            % For dc apparatuses: smaller
            Pos_Apparatus{i} = Pos_Bus{Bus} + Shift_Apparatus;
            set_param(FullName_Apparatus{i},'position',[Pos_Apparatus{i},Pos_Apparatus{i}+Size_Apparatus-[0,20]]);
        elseif 2000<=ApparatusType{i} && ApparatusType{i}<=2090
            % For interlink apparatuses: larger
            Pos_Apparatus{i} = Pos_Bus{Bus(1)} + Shift_Apparatus;
            set_param(FullName_Apparatus{i},'position',[Pos_Apparatus{i},Pos_Apparatus{i}+Size_Apparatus+[0,40]]);
        end
        set_param(FullName_Apparatus{i},'Orientation','left');
        
        % Set common variables
      	set_param(gcb,'Sbase','Sbase');
        set_param(gcb,'Vbase','Vbase');
        set_param(gcb,'Ts','Ts');
        
        % For ac apparatus only
        if (ApparatusType{i}<1000) || (2000<=ApparatusType{i} && ApparatusType{i}<=2090)
            set_param(gcb,'Wbase','Wbase');
        end
        
        % For active apparatus only
        if (0<=ApparatusType{i} && ApparatusType{i}<90) || ...
           (1000<=ApparatusType{i} && ApparatusType{i}<1090) || ...
           (2000<=ApparatusType{i} && ApparatusType{i}<2090)
            
            % Set system object parameters
            set_param(gcb,'ApparatusType',['ApparatusType{' num2str(i) '}']);
            set_param(gcb,'ApparatusPara',['ApparatusPara{' num2str(i) '}']);
            set_param(gcb,'PowerFlow',['ApparatusPowerFlow{' num2str(i) '}']);
            set_param(gcb,'x0',['x_e{' num2str(i) '}']);
            set_param(gcb,'OtherInputs',['OtherInputs{' num2str(i) '}']);

            % Set discretization method
            switch DiscreMethod
                case 1
                    ApparatusDiscreMethod = 'Forward Euler';
                case 2
                    ApparatusDiscreMethod = 'Hybrid Trapezoidal';
                case 3
                    ApparatusDiscreMethod = 'Virtual Damping';
                otherwise
                    error(['Error: Wrong discretization method.'])
            end
            set_param(gcb,'DiscreMethod',ApparatusDiscreMethod);
            set_param(gcb,'LinearizationTimes',num2str(LinearizationTimes));
            if DirectFeedthrough == 1
                set_param(gcb,'DirectFeedthrough','on');
            else
                set_param(gcb,'DirectFeedthrough','off');
            end
            set_param(gcb,'EnableInsideModification','on');
            if DiscreDampingFlag == 1
                set_param(gcb,'DiscreDampingFlag','on');
                set_param(gcb,'DiscreDampingValue',['ApparatusDiscreDamping{' num2str(i) '}']);
            else
                set_param(gcb,'DiscreDampingFlag','off');
            end
            set_param(gcb,'EnableInsideModification','off');
            
        end
        
        % If the apparatus is an infinite bus
        if ApparatusType{i} == 0090        % Ac infinite bus
            set_param(gcb,'vd',['PowerFlow{' num2str(i) '}(3)']);
            set_param(gcb,'theta0',['PowerFlow{' num2str(i) '}(4)']);
            set_param(gcb,'w','Wbase');
        elseif ApparatusType{i} == 1090    % Dc infinite bus
            set_param(gcb,'v',['PowerFlow{' num2str(i) '}(3)']);
        end
        
    end
end

end