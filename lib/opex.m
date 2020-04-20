classdef opex
    % OPEX Library for all OPEX-related equations
    
    %% PROPERTIES
    properties
        ...
    end

    %% CONSTANT PROPERTIES
    properties (Constant)
        CEPCI19 = 607.5; % CEPCI index for 2019 (annual)
        CEPCI18 = 603.1; % CEPCI index for 2018 (annual)
        CEPCI01 = 394.3; % CEPCI index for 2001 (annual)
        
        USSG19 = 1.3493; % USD to SGD forex rate for 2019 (annual)
        USSG18 = 1.3912; % USD to SGD forex rate for 2018 (annual)
        USSG01 = 1.7912; % USD to SGD forex rate for 2001 (annual)
        
        CPISG19 = 115.0; % SG consumer price index for 2019 (annual, SG benchmark 2010 = 100)
        CPISG18 = 113.8; % SG consumer price index for 2018 (annual, SG benchmark 2010 = 100)
        CPISG01 = 86.05; % SG consumer price index for 2001 (annual, SG benchmark 2010 = 100)
        CPIUS19 = 257.0; % US consumer price index for 2019 (annual, US benchmark 1983 = 100)
        CPIUS18 = 251.2; % US consumer price index for 2018 (annual, US benchmark 1983 = 100)
        CPIUS01 = 176.7; % US consumer price index for 2001 (annual, US benchmark 1983 = 100)
        
        runtime = 8000; % Operational runtime (h/yr)
        shiftdur = 8; % Duration per workshift (h/shift)
        shiftperwk = 5; % Number of workshifts per week
        yearww = 49; % Number of work weeks per year
        SF = 8000/(365*24); % Stream factor
        
        % Utility Costs (per GJ basis)
        % From "Analysis, Synthesis & Design of Chemical Processes, 5th Ed. by Turton et al."
        utilLPS = 4.54; % Utility cost for LPS (5 barg, 160 degC) ($/GJ)
        utilMPS = 4.77; % Utility cost for MPS (10 barg, 184 degC) ($/GJ)
        utilHPS = 5.66; % Utility cost for HPS (41 barg, 254 degC) ($/GJ)
        utilCW = 0.378; % Utility cost for cooling water (30-45 degC) ($/GJ)
        utilChW = 4.77; % Utility cost for chilled water (5 degC) ($/GJ)
        utilLTR = 8.49; % Utility cost for low temperature refrigerant (-20 degC) ($/GJ)
        utilVLTR = 14.12; % Utility cost for very low temperature refrigerant (-50 degC) ($/GJ)
        utilelec = 18.72; % Utility cost for electricity (110-440 V) ($/GJ)
    end
    
    %% METHODS
    methods
        ...
    end

    %% STATIC METHODS
    methods (Static)
        function NOL = operatorspershift(P,Nnp)
            % OPERATORSPERSHIFT Calculate number of operators per shift
            %   NOL = OPERATORSPERSHIFT(P,Nnp) calculates the number of
            %   operators per shift, where:
            %   - P = number of processing steps involving particulate
            %   solids (P=0 for fluid-processing plants)
            %   - Nnp = number of non-particulate/fluid handling
            %   equipment/steps (include compressors, towers, reactors,
            %   heaters and exchangers; exclude pumps, vessels and tanks)
            %
            %   Returns:
            %   - NOL = number of operators required per shift
            
            NOL = round(sqrt(6.29 + 31.7*P^2 + 0.23*Nnp));
        end
        
        function COL = labourcost(NOL,wage)
            % LABOURCOST Calculate annualised labour cost
            % COL = LABOURCOST(NOL,wage) calculates the annualised labour
            % cost, where:
            %   - NOL = total number of operators required
            %   - wage = annualised per-operator wage ($)
            %
            %   Returns:
            %   - COL = annualised labour cost ($)
            
            shiftperyr = opex.runtime / opex.shiftdur;
            shiftperopperyr = opex.yearww * opex.shiftperwk;
            Nop = round(shiftperyr / shiftperopperyr * NOL);
            COL = Nop * wage;
        end
        
        function CUT = costofutil(HPS,MPS,LPS,CW,ChW,LTR,VLTR,elec)
            % COSTOFUTIL Calculate cost of utilities
            % CUT = COSTOFUTIL(utilities) calculates the annualised cost of
            % utilities, where:
            %   - HPS = annual consumption of high-pressure steam (GJ)
            %   - MPS = annual consumption of medium-pressure steam (GJ)
            %   - LPS = annual consumption of low-pressure steam (GJ)
            %   - CW = annual consumption of cooling water (GJ)
            %   - ChW = annual consumption of chilled water (GJ)
            %   - LTR = annual consumption of low-temperature refrigerant (GJ)
            %   - VLTR = annual consumption of very low-temperature refrigerant (GJ)
            %   - elec = annual consumption of electricity (GJ)
            %
            %   Returns:
            %   - CUT = annualised cost of utilities ($)
            
            CUT = [opex.utilHPS,opex.utilMPS,opex.utilLPS,...
                opex.utilCW,opex.utilChW,...
                opex.utilLTR,opex.utilVLTR,opex.utilelec]*...
                [HPS,MPS,LPS,CW,ChW,LTR,VLTR,elec]';
        end
        
        function [COMd,COM,d,DMC,FMC,GE] = costofmanfc(FCI,COL,CRM,CWT,CUT)
            % COSTOFMANFC Calculate all components of annualised total cost of manufacture (COM)
            % [COMd,COM,d,DMC,FMC,GE] = COSTOFMANFC(FCI,COL,CRM,CWT,CUT)
            % calculates all components of annualised cost of manufacture
            % (COM), where:
            %   - FCI = fixed capital investment (equals grassroots cost
            %   CGR for greenfield projects, or total module cost CTM for
            %   brownfield projects) ($)
            %   - COL = annualised cost of operating labour ($)
            %   - CRM = annualised cost of raw materials ($)
            %   - CWT = annualised cost of waste treatment ($)
            %   - CUT = annualised cost of utilities ($)
            %
            %   Returns:
            %   - COMd = total annualised cost of manufacturing, without
            %   depreciation ($)
            %   - COM = total annualised cost of manufacturing, with
            %   depreciation ($)
            %   - d = annualised depreciation, assumed to approximate with
            %   10% of FCI ($)
            %   - DMC = annualised direct manufacturing cost, which is
            %   dependent on production rate ($)
            %   - FMC = annualised fixed manufacturing cost, which is not
            %   dependent on production rate ($)
            %   - GE = annualised general manufacturing expense, which are
            %   costs associated with management level and administrative
            %   activities not directly related to the manufacturing process
            
            COMd = 0.18 * FCI + 2.73 * COL + 1.23 * (CRM + CWT + CUT);
            d = 0.1 * FCI;
            COM = COMd + d;
            DMC = CRM + CWT + CUT + 1.33 * COL + 0.069 * FCI + 0.03 * COM;
            FMC = 0.708 * COL + 0.068 * FCI;
            GE = 0.177 * COL + 0.009 * FCI + 0.16 * COM;
        end
    end
end