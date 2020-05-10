classdef capex
    % CAPEX Library for all CAPEX-related equations
    
    %% PROPERTIES
    properties
        ...
    end

    %% CONSTANT PROPERTIES
    properties (Constant)
        
        % CEPCI Index for 20yy (annual)
        % To access, e.g. capex.CEPCI(2019)
        CEPCI = containers.Map([2001, 2018, 2019], [394.3, 603.1, 607.5])
        
        % USD to SGD forex rate for 20yy (annualised average)
        % To access, e.g. capex.USSG(2019)
        USSG = containers.Map([2001, 2018, 2019], [1.7912, 1.3912, 1.3493])
        
        % Consumer price index for 20yy (annual)
        % To access, e.g. CPI('SG')(2019)
        CPI = containers.Map(['SG', 'US'], ...
            [containers.Map([2001, 2016, 2018, 2019], [86.05, 112.6, 113.8, 115.0]), ...
            containers.Map([2001, 2016, 2018, 2019], [176.7, 240.0, 251.2, 257.0])])        
        
    end
    
    %% METHODS
    methods
        ...
    end

    %% STATIC METHODS
    methods (Static)
        function Cpo = eqptpurcost(A,K1,K2,K3)
            % EQPTPURCOST Calculate equipment purchased cost (Cp^o)
            %   Cp = EQPTPURCOST(A,K1,K2,K3) calculates equipment purchased
            %   cost at ambient pressure and using carbon steel as MOC,
            %   where:
            %   - A = equipment capacity (various units)
            %   - K1, K2, K3 = cost correlation factors
            %
            %   Returns:
            %   - Cpo = equipment purchased cost ($)
            
            Cpo = 10 .^ (K1 + K2.*log10(A) + K3.*(log10(A)).^2);
        end
        
        function FP = pressurefacves(D,ts,P)
            % PRESSUREFACVES Calculate pressure factor (F_P) for vessels
            %   CBM = PRESSUREFACVES(D,P) calculates bare module factor for
            %   vessels at specified elevated pressure and MOC,
            %   where:
            %   - D = vessel diameter (m)
            %   - ts = vessel thickness (in)
            %   - P = pressure (barg)
            %
            %   Returns:
            %   - FP = amplification factor for pressure
            if P < -0.5
                FP = 1.25;
            elseif P > -0.5 && ts < dsg.tmin
                FP = 1;
            else
                FP = max(((P+1)*D / (2*(850-0.6*(P+1))) + 0.00315) / 0.0063, 1);
            end
        end
        
        function FP = pressurefacanc(P,C1,C2,C3)
            % PRESSUREFACANC Calculate pressure factor (F_P) for anciliary equipment (e.g. pumps and exchangers)
            %   CBM = PRESSUREFACANC(P,C1,C2,C3) calculates bare module factor
            %   for anciliary equipment (e.g. pumps and exchangers) at
            %   specified elevated pressure and MOC, where:
            %   - P = pressure (barg)
            %   - C1, C2, C3 = pressure correlation factors
            %
            %   Returns:
            %   - FP = amplification factor for pressure
            
            FP = capex.eqptpurcost(P,C1,C2,C3);
            % borrowing quadratic-exponential relation
        end
        
        function FBM = baremodfac(B1,B2,FM,FP)
            % BAREMODFAC Calculate bare module factor (F_BM)
            %   CBM = BAREMODFAC(B1,B2,FM,FP) calculates
            %   bare module factor at specified elevated pressure and MOC, where:
            %   - B1, B2 = bare module correlation factors
            %   - FM = amplification factor for material of construction
            %   (MOC)
            %   - FP = amplification factor for pressure
            %
            %   Returns:
            %   - FBM = bare module factor (dimensionless)
            
            FBM = B1 + B2.*FM.*FP;
        end
        
        function CBM = baremodcost(Cpo, FBM)
            % BAREMODCOST Calculate bare module cost (C_BM)
            %   CBM = BAREMODCOST(Cpo, FBM) calculates
            %   bare module cost at specified elevated pressure and MOC,
            %   where:
            %   - Cpo = equipment purchased cost ($)
            %   - FBM = bare module factor (dimensionless)
            %
            %   Returns:
            %   - CBM = bare module cost ($)
            
            CBM = FBM .* Cpo;
        end
        
        function CTM = totmodcost(CBM)
            % TOTMODCOST Calculate total module cost (C_TM)
            %   CBM = TOTMODCOST(Cpo, FBM) calculates
            %   total module cost at specified elevated pressure and MOC,
            %   where:
            %   - CBM = bare module cost ($)
            %
            %   Returns:
            %   - CTM = total module cost ($)
            
            CTM = 1.18 .* CBM;
        end
        
        function CGR = grasscost(CTM,Cpo)
            % GRASSCOST Calculate grassroots cost (CGR)
            %   CGR = GRASSCOST(CTM,Cpo) calculates grasscosts cost
            %   at specified elevated pressure and MOC,
            %   where:
            %   - CTM = total module cost ($)
            %   - Cpo = purchased equipment cost at ambient pressure and
            %   carbon steel MOC ($)
            %
            %   Returns:
            %   - CGR = grassroots cost ($)
            
            CGR = CTM + 0.5 .* Cpo;
        end
    end
end