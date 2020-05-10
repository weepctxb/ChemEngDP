classdef dsg
    % DSG Library for all mechanical design-related equations
    
    %% PROPERTIES
    properties
        ...
    end

    %% CONSTANT PROPERTIES
    properties (Constant)
        Patm = 14.696 % Standard atmospheric pressure (psi)
        Patmb = 1.01325 % Standard atmospheric pressure (bar)
        Troom = 77 % Ambient temperature (degF)
        tmin = 1/4 % Universal minimum allowable vessel thickness (in)
        tc = 0.125 % Corrosion allowance (in) for both corrosive and non-corrosive conditions (default is 1/8)
        rho = 0.2836 % Density of SA-285C/SA-387B/carbon/low-alloy steels (lb/in^3)
        
        g = 9.80665; % standard Earth gravitational acceleration (m/s^2)
        R = 8.31446261815324; % universal ideal gas constant (J/(K.mol))
        Ta = 10; % minimum heat exchanger temperature approach (K)
    end
    
    %% METHODS
    methods
        ...
    end

    %% STATIC METHODS
    methods (Static)
        function Pd = designP(Po)
            % DESIGNP Calculate design pressure for pressure vessels
            %   Pd = DESIGNP(Po) calculates design pressure based on most
            %   deviated operating pressure for both pressure and vacuum
            %   vessels, where:
            %   - Po = most deviated operating pressure (psig)
            %
            %   Returns:
            %   - Pd = design pressure (psig)
            
            if 0<=Po && Po<=5
                Pd = 10;
            elseif 5<Po && Po<10
                Pd = max(10,...
                    exp(0.60608+0.91615*log(Po)+0.0015655*(log(Po))^2));
            elseif 10<=Po && Po<=1e3
                Pd = exp(0.60608+0.91615*log(Po)+0.0015655*(log(Po))^2);
            elseif 1e3<Po
                Pd = 1.1*Po;
            elseif Po<0
                Pd = Po;
            else
                error("Input Po out of supported range!")
            end
        end
        
        function Td = designT(To,varargin)
            % DESIGNT Calculate design temperature for pressure and vacuum vessels
            %   Td = DESIGNT(To, heuristic) calculates design temperature
            %   based on most deviated operating temperature for both cases
            %   above and below ambient temperature, where:
            %   - To = most deviated operating temperature (degF)
            %   - heuristic = (either "Towler" or "Turton", optional
            %   string)
            %
            %   Caution: If heuristic is not specified, Towler's design
            %   temperature heuristic will be used by default.
            %
            %   Returns:
            %   - Td = design pressure (degF)
            
            if nargin == 1
                heuristic = "towler";
            else
                varargin{1}
                heuristic = lower(varargin{1});
            end
            
            if contains(heuristic,"towler")
                if To<dsg.Troom
                    Td = To - 25;
                else
                    Td = To + 50;
                end
            elseif contains(heuristic,"turton")
                if -22<=To && To<=644
                    Td = To + 45;
                else
                   error("To temperature input out of supported " + ...
                       "range using Turton heuristic!")
                end
            else
                error("Heuristic not supported! " + ...
                    "Please check heuristic input!")
            end
        end
        
        function [Smax, MOC] = maxstress(Td, varargin)
            % MAXSTRESS Calculate maximum allowable stress for pressure vessel material (for SA-317L)
            %   [Smax, MOC] = MAXSTRESS(Td, MOC) calculates maximum allowable
            %   stress for pressure vessels based on design temperature, where:
            %   - Td = design temperature (degF)
            %   - MOC = user-specified material of construction (optional input)
            %
            %   Returns:
            %   - Smax = maximum allowable stress for pressure vessel MOC (psi)
            %   - MOC = prescribed material of construction which is
            %   in stainless steel family (string, optional). If MOC is not
            %   user-specified, the returned MOC will be a default value
            %   (SA-285C or SA-387B).
            %
            %   See also: elasmod (vacuum vessels)
            
            if nargin > 1
                MOC = varargin{1};
            else
                MOC = "";
            end
            
            if contains(MOC,"317L")
                if -20<=Td & Td<=68
                    MOC = "SA-317L";
                    Smax = 25286;
                elseif 68<Td & Td<200
                    MOC = "SA-317L";
                    Smax = 22957;
                elseif 200<=Td & Td<400
                    MOC = "SA-317L";
                    Smax = 20957;
                elseif 400<=Td & Td<600
                    MOC = "SA-317L";
                    Smax = 19400;
                elseif 600<=Td & Td<=800
                    MOC = "SA-317L";
                    Smax = 17633;
                elseif 800<=Td & Td<=1000
                    MOC = "SA-317L";
                    Smax = 16733;
                elseif 1000<=Td & Td<=1200
                    MOC = "SA-317L";
                    Smax = 15767;
                elseif 1200<=Td & Td<=1400
                    MOC = "SA-317L";
                    Smax = 12857;
                elseif 1400<=Td & Td<=1600
                    MOC = "SA-317L";
                    Smax = 8300;
                else
                    error("Input Td out of supported range!")
                end
            else
                if -20<=Td & Td<=650
                    MOC = "SA-285C";
                    Smax = 13750;
                elseif 650<Td & Td<750
                    MOC = "SA-387B"
                    Smax = 15000;
                elseif 750<=Td & Td<800
                    MOC = "SA-387B"
                    Smax = 14750;
                elseif 800<=Td & Td<850
                    MOC = "SA-387B"
                    Smax = 14200;
                elseif 850<=Td & Td<=900
                    MOC = "SA-387B"
                    Smax = 13100;
                else
                    error("Input Td out of supported range!")
                end
            end
        end
        
        function [EM, MOC] = elasmod(Td,MOC)
            % ELASMOD Calculate modulus of elasticity for vacuum vessel material
            %   EM = ELASMOD(Td,MOC) calculates the modulus of elasticity
            %   for vacuum vessels based on design temperature and
            %   specified material of construction (MOC), where:
            %   - Td = design pressure (degF)
            %   - MOC = material of construction (only either "carbon" or
            %   "low-alloy", string)
            %
            %   Returns:
            %   - EM = modulus of elasticity for vacuum vessel MOC (psi)
            %   - MOC = returns the input MOC for consistency with the
            %   equivalent computation for pressure vessels (optional string)
            %
            %   See also: MAXSTRESS (pressure vessels)
            
            if contains(lower(MOC),"carbon")
                MOC = "carbon";
                if Td<=-20
                    EM = 30.2e6;
                elseif -20<Td && Td<=200
                    EM = 29.5e6;
                elseif 200<Td && Td<=400
                    EM = 28.3e6;
                elseif 400<Td && Td<=650
                    EM = 26.0e6;
                else
                    error("Td out of supported range for carbon steel!")
                end
            elseif contains(lower(MOC),"low") && ...
                    contains(lower(MOC),"alloy")
                MOC = "low-alloy";
                if Td<=-20
                    EM = 30.2e6;
                elseif -20<=Td && Td<200
                    EM = 29.5e6;
                elseif 200<Td && Td<=400
                    EM = 28.6e6;
                elseif 400<Td && Td<=650
                    EM = 27.0e6;
                elseif 650<Td && Td<=700
                    EM = 26.6e6;
                elseif 700<Td && Td<=800
                    EM = 25.7e6;
                elseif 800<Td && Td<=900
                    EM = 24.5e6;
                else
                    error("Td out of supported range for low-alloy steel!")
                end
            else
                error("Specified MOC not found! Please check MOC input!")
            end
        end
        
        function [tp, E] = wallthk(Pd,Di,Smax)
            % WALLTHK Calculate cylindrical shell wall thickness for pressure vessels
            %   [tp, E] = WALLTHK(Pd,Di,Smax) calculates cylindrical shell
            %   wall thickness for pressure vessels, including minimum
            %   thickness check for structural rigidity, where:
            %   - Pd = design pressure (psig)
            %   - Di = internal diameter (in)
            %   - Smax = maximum allowable stress (psi)
            %
            %   Returns:
            %   - tp = cylindrical shell wall thickness for pressure
            %   vessels (in)
            %   - E = fractional weld efficiency used (string, optional)
            
            E = 0.85; % First assume 10% X-ray spot check
            tp = Pd*Di / (2*Smax*E - 1.2*Pd);
            
            if tp>1.25 % If tp not large enough, 100% X-ray check needed
                E = 1;
                tp = Pd*Di / (2*Smax*E - 1.2*Pd);
            end
            
            % Check if minimum wall thickness to provide rigidity satisfied
            if Di<48
                tp = max(1/4, tp);
            elseif 48<=Di && Di<72
                tp = max(5/16, tp);
            elseif 72<=Di && Di<96
                tp = max(3/8, tp);
            elseif 96<=Di && Di<120
                tp = max(7/16, tp);
            elseif 120<=Di && Di<144
                tp = max(1/2, tp);
            else
                warning("Input Di out of supported range for" + ...
                    " minimum tp check. Assuming minimum tp = 1/2 in!")
                tp = max(1/2, tp);
            end
        end
        
        function [tp, tE, tEC] = wallthkvac(Pd,Do,Di,L,EM)
            % WALLTHKVAC Calculate cylindrical shell wall thickness for vacuum vessels
            %   [tp, tE, tEC] = WALLTHKVAC(Pd,Do,Di,L,EM) calculates
            %   cylindrical shell wall thickness for vacuum vessels, where:
            %   - Pd = design pressure (psig)
            %   - Do = external diameter (in)
            %   - Di = internal diameter (in)
            %   - L = vessel length (in)
            %   - EM = modulus of elasticity (psi)
            %
            %   Caution: This assumes that the value of Do, which is
            %   dependent on tp itself, is known.
            %
            %   Returns:
            %   - tp = cylindrical shell wall thickness for vacuum vessels 
            %   (in)
            %   - tE = necessary thickness for vacuum vessels (in, optional)
            %   - tEC = correction factor for vacuum vessels (in, optional)
            
            tE = 1.3 * Do * (Pd*L / (EM*Do)) ^ 0.4;
            if tE/Do > 0.05
                warning("tE is > 0.05*Do, which is" + ...
                    " outside the validity range for tE computation!" + ...
                    " Nevertheless carrying on with calculation - beware!")
            end
            tEC = L * (0.18*Di - 2.2)*1e-5 - 0.19;
            tp = tE + tEC;
        end
        
        function ts = shellthkhorz(tp)
            % SHELLTHKHORZ Calculate shell thickness for horizontal vessels
            %   ts = SHELLTHKHORZ(tp) calculates the final shell thickness
            %   for horizontal vessels after incorporating corrosion
            %   allowance, where:
            %   - tp = wall thickness (in)
            %
            %   Returns:
            %   - ts = shell thickness with corrosion allowance for
            %   horizontal vessels (in)
            
            ts = tp + dsg.tc; % add corrosion allowance
        end
        
        function tw = windalw(Do,L,Smax)
            % WINDALW Calculate wind/earthquake allowance for vertical vessels
            %   tw = WINDALW(Do,L,Smax) calculates the additional thickness
            %   allowance for wind/earthquake (only for vertical
            %   thickness!), where:
            %   - Do = external diameter (in)
            %   - L = internal tangent-to-tangent height (in)
            %   - Smax = maximum allowable stress (psi)
            %
            %   Returns:
            %   - tw = wind/earthquake allowance for vertical vessels (in)
            %
            %   Caution: Using WINDALW requires an assumed value of Do
            %   which is dependent on tw. If Do is unknown, use
            %   SHELLTHKVERT directly instead which internally calls
            %   WINDALW.
            %   
            %   See also: SHELLTHKVERT
            
            tw = 0.22 * (Do+18) * L^2 / (Smax*Do^2);
        end
        
        function [ts, tv, tw] = shellthkvert(tp,Di,L,Smax)
            % SHELLTHKVERT Calculate shell thickness for vertical vessels
            %   [ts, tv] = SHELLTHKVERT(tp,Di,L,Smax) calculates the final
            %   shell thickness for vertical vessels after incorporating
            %   corrosion allowance and wind/earthquake allowance, where:
            %   - tp = wall thickness (in)
            %   - Di = internal diameter (in)
            %   - L = internal tangent-to-tangent height (in)
            %   - Smax = maximum allowable stress of MOC (psi)
            %
            %   Caution: For vertical vessels, ts is implicit, i.e.
            %   ts=f(Do(ts)). This implementation also calls WINDALW to
            %   calculate the wind/earthquake allowance for vertical
            %   vessels, and the computation is iterated until convergence
            %   in ts is reached.
            %
            %   Returns:
            %   - ts = shell thickness with wind allowance after adding
            %   corrosion allowance for vertical vessels (in)
            %   - tv = shell thickness with wind allowance before adding
            %   corrosion allowance for vertical vessels (in, optional)
            %   - tw = wind allowance (in, optional)
            %
            %   See also: WINDALW
            
            ts0 = 2*tp; % dummy initialisation
            Do = Di + 2*ts0;
            tw = dsg.windalw(Do,L,Smax);
            tv = (tp + (tp+tw)) / 2;
            ts1 = tv + dsg.tc; % add corrosion allowance
            
            reltol = 1e-9; i = 0;
            while abs(ts1-ts0)/ts0 > reltol && i < 1e3
                ts0 = ts1; i = i + 1;
                Do = Di + 2*ts0;
                tw = dsg.windalw(Do,L,Smax);
                tv = (tp + (tp+tw)) / 2;
                ts1 = tv + dsg.tc; % add corrosion allowance
            end
            
            if i == 1e3
                warning("Vertical vessel thickness failed to converge!" + ...
                    " Nevertheless carrying on with calculation - beware!")
            end
            
            ts = ts1;
        end
        
        function tsfinal = ceilplatethk(ts)
            % CEILPLATETHK Round up metal plate thickness to nearest increment
            %   tsfinal = CEILPLATETHK(ts) rounds up shell wall thickness
            %   ts to the nearest industrial increment for metal plates,
            %   where:
            %   - ts = shell wall thickness before before rounding to
            %   nearest increment (in)
            %
            %   Returns:
            %   - tsfinal = final shell wall thickness after rounding to
            %   nearest increment (in)
            
            if 3/16<=ts && ts<=1/2
                acc = 1/16;
            elseif 1/2<ts && ts<=2
                acc = 1/8;
            elseif 2<ts && ts<=3
                acc = 1/4;
            elseif ts<3/16
                error("ts too low (<3/16 in). Please debug calculations!")
            else
                warning("Calculate ts not in supported range." + ...
                    " Assuming metal plate thickness in increments"+ ...
                    " of 1/4 inches above 3 inches!")
                acc = 1/4;
            end
            
            tsfinal = ceil(ts/acc)*acc;
        end
        
        function W = vesselweight(Di,tsfinal,L,rho)
            % VESSELWEIGHT Calculate final weight of vessel
            %   W = VESSELWEIGHT(Di,ts,L,rho) calculates the final weight
            %   of the vessel with the shell and two 2:1 elliptical heads,
            %   where:
            %   - Di = internal diameter (in)
            %   - tsfinal = shell thickness with corrosion allowance,
            %   rounded to nearest thickness increment for metal plates (in)
            %   - L = internal tangent-to-tangent length/height (in)
            %   - rho = density of material of construction (MOC) (lb/in^3)
            %
            %   For cost estimation purposes the head thickness is assumed
            %   to be equal to shell thickness.
            %
            %   Returns:
            %   - W = weight of vessel (lb)
            
            W = pi * (Di+tsfinal) * (L+0.8*Di) * tsfinal * rho;
        end
        
        function V = vesselvol(Do,L)
            % VESSELVOL Calculate final weight of vessel
            %   V = VESSELVOL(Do,L) calculates the final volume
            %   of the vessel with the shell and two 2:1 elliptical heads,
            %   where:
            %   - Do = external diameter (in)
            %   - L = internal tangent-to-tangent length/height (in)
            %
            %   For cost estimation purposes the head thickness is assumed
            %   to be equal to shell thickness.
            %
            %   Returns:
            %   - V = volume of vessel (in^3)
            
            Vcyl = pi * Do^2 / 4 * L;
            H = Do / 4;
            Vheads = 4/3 * pi * H * (Do/2)^2;
            V = Vcyl + Vheads;
        end
        
        function [md, in] = designhorzpres(Po,To,Di,L)
            % DESIGNHORZPRES The main function to be called for designing horizontal pressure vessels
            %   mechdesign = DESIGNHORZPRES(Po,To,Di,L) is the main
            %   function to be called for designing horizontal pressure
            %   vessels, where:
            %   - Po = most deviated operating pressure from ambient
            %   pressure (psig)
            %   - To = most deviated operating temperature from ambient
            %   temperature (degF)
            %   - Di = internal diameter (in)
            %   - L = tangent-to-tangent horizontal length (in)
            %
            %   Returns and displays a data struct "mechdesign" consisting
            %   of:
            %   - Pd = design pressure (psig)
            %   - Td = design pressure (degF)
            %   - MOC = material of construction to use (string)
            %   - Smax = maximum allowable stress of MOC used (psi)
            %   - E = weld efficiency to use (dimensionless)
            %   - tp = wall thickness (in)
            %   - tc = corrosion allowance used (= 1/8 in)
            %   - ts = tp with tc (in)
            %   - tsfinal = ts rounded up to next increment in metal plate thickness (in)
            %   - Do = external diameter (in)
            %   - W = total vessel weight (lb)
            %   - V = total vessel volume (in^3)
            %
            %   Example implementation:
            %   mechdesign1 = dsg.designhorzpres(470,850,78,480)
            
            in = []; input = []; md = []; mechdesign = [];
            
            md.Pd = dsg.designP(Po);
            md.Td = dsg.designT(To);
            [md.Smax, md.MOC] = dsg.maxstress(md.Td);
            [md.tp, md.E] = dsg.wallthk(md.Pd, Di, md.Smax);
            md.tc = dsg.tc;
            md.ts = dsg.shellthkhorz(md.tp);
            md.tsfinal = dsg.ceilplatethk(md.ts);
            md.Do = Di + 2 * md.tsfinal;
            md.W = dsg.vesselweight(Di, md.tsfinal, L, dsg.rho);
            md.V = dsg.vesselvol(md.Do, L);
            
            in.Po = Po;
            in.To = To;
            in.Di = Di;
            in.L = L;
        end
        
        function [md, in] = designvertpres(Po,To,Di,L)
            % DESIGNVERTPRES The main function to be called for designing vertical pressure vessels
            %   mechdesign = DESIGNVERTPRES(Po,To,Di,L) is the main
            %   function to be called for designing vertical pressure
            %   vessels, where:
            %   - Po = most deviated operating pressure from ambient
            %   pressure (psig)
            %   - To = most deviated operating temperature from ambient
            %   temperature (degF)
            %   - Di = internal diameter (in)
            %   - L = tangent-to-tangent horizontal height (in)
            %
            %   Returns and displays a data struct "mechdesign" consisting
            %   of:
            %   - Pd = design pressure (psig)
            %   - Td = design pressure (degF)
            %   - MOC = material of construction to use (string)
            %   - Smax = maximum allowable stress of MOC used (psi)
            %   - E = weld efficiency to use (dimensionless)
            %   - tp = wall thickness (in)
            %   - tc = corrosion allowance used (= 1/8 in)
            %   - tw = wind/earthquake allowance for vertical vessels (in)
            %   - tv = tp with tw without tc (in)
            %   - ts = tp with both tw and tc (in)
            %   - tsfinal = ts rounded up to next increment in metal plate thickness (in)
            %   - Do = external diameter (in)
            %   - W = total vessel weight (lb)
            %   - V = total vessel volume (in^3)
            %
            %   Example implementation:
            %   mechdesign2 = dsg.designvertpres(95.5,150,120,2544)
            
            in = []; input = []; md = []; mechdesign = [];
            
            md.Pd = dsg.designP(Po);
            md.Td = dsg.designT(To);
            [md.Smax, md.MOC] = dsg.maxstress(md.Td);
            [md.tp, md.E] = dsg.wallthk(md.Pd, Di, md.Smax);
            md.tc = dsg.tc;
            [md.ts, md.tv, md.tw] = dsg.shellthkvert(md.tp, Di, L, md.Smax);
            md.tsfinal = dsg.ceilplatethk(md.ts);
            md.Do = Di + 2 * md.tsfinal;
            md.W = dsg.vesselweight(Di, md.tsfinal, L, dsg.rho);
            md.V = dsg.vesselvol(md.Do, L);
            
            in.Po = Po;
            in.To = To;
            in.Di = Di;
            in.L = L;
        end
        
        function [md, in] = designvac(Po,To,Di,L)
            % DESIGNVAC The main function to be called for designing vacuum vessels
            %   mechdesign = DESIGNVAC(Po,To,Di,L) is the main
            %   function to be called for designing vacuum vessels, where:
            %   - Po = most deviated operating pressure from ambient
            %   pressure (psig)
            %   - To = most deviated operating temperature from ambient
            %   temperature (degF)
            %   - Di = internal diameter (in)
            %   - L = tangent-to-tangent horizontal length/height (in)
            %
            %   Returns and displays a data struct "mechdesign" consisting
            %   of:
            %   - Pd = design pressure (psig)
            %   - Td = design pressure (degF)
            %   - MOC = material of construction to use (string)
            %   - EM = modulus of elasticity of MOC used (psi)
            %   - tE = vacuum wall thickness (in)
            %   - tEC = vacuum wall correction factor (in)
            %   - tp = tE with tEC (in)
            %   - tc = corrosion allowance used (= 1/8 in)
            %   - ts = tp with tc (in)
            %   - tsfinal = ts rounded up to next increment in metal plate thickness (in)
            %   - Do = external diameter (in)
            %   - W = total vessel weight (lb)
            %   - V = total vessel volume (in^3)
            %
            %   Example implementation:
            %   mechdesign3 = dsg.designvac(7.977,257,168,1080)
            
            in = []; input = []; md = []; mechdesign = [];
            
            md.Pd = dsg.Patm - Po;
            md.Td = dsg.designT(To);
            if -20<md.Td && md.Td<=650
                [md.EM, md.MOC] = dsg.elasmod(md.Td,"carbon");
            elseif 650<md.Td && md.Td<900
                [md.EM, md.MOC] = dsg.elasmod(md.Td,"low-alloy");
            else
                error("Td out of supported range for both" + ...
                    " carbon and low-alloy steel!")
            end
            
            ts0 = 1; % dummy initialisation
            md.Do = Di + 2*ts0;
            [md.tp, md.tE, md.tEC] = dsg.wallthkvac(md.Pd, md.Do, Di, L, md.EM);
            md.tc = dsg.tc;
            ts1 = dsg.shellthkhorz(md.tp); % horz/vert orientation does not matter for vacuum
            
            reltol = 1e-9; i = 0;
            while abs(ts1-ts0)/ts0 > reltol && i < 1e3
                ts0 = ts1; i = i + 1;
                md.Do = Di + 2*ts0;
                [md.tp, md.tE, md.tEC] = dsg.wallthkvac(md.Pd, md.Do, Di, L, md.EM);
                md.tc = dsg.tc;
                ts1 = dsg.shellthkhorz(md.tp);
            end
            
            if i == 1e3
                warning("Vacuum vessel thickness failed to converge!" + ...
                    " Nevertheless carrying on with calculation - beware!")
            end
            
            md.ts = ts1;
            
            md.tsfinal = dsg.ceilplatethk(md.ts);
            md.Do = Di + 2 * md.tsfinal;
            md.W = dsg.vesselweight(Di, md.tsfinal, L, dsg.rho);
            md.V = dsg.vesselvol(md.Do, L);
            
            in.Po = Po;
            in.To = To;
            in.Di = Di;
            in.L = L;
        end
        
        function [comppower, compeff, T2] = sizecompressor(m,P1,P2,T1,cp,cv,varargin)
            % SIZECOMPRESSOR Conducts compressor sizing by determining required compressor power
            %   [comppower, compeff, T2] = sizecompressor(m,P1,P2,T1,cp,cv,Z)
            %   calculates poewr required to operate a compressor based on
            %   its flow rate and inlet/outlet pressures, where:
            %   - m = mass flow rate through compressor (kg/h)
            %   - P1 = gas inlet pressure (any pressure units)
            %   - P2 = gas inlet pressure (same pressure unit as P1)
            %   - T1 = gas inlet temperature (K)
            %   - cp = constant-pressure heat capacity of gas
            %   - cv = constant-volume heat capacity of gas
            %   - Z = gas compressibility factor (optional input, if not
            %   specified the default value of one is assumed)
            %
            %   Returns:
            %   - comppower = required compressor power (kW)
            %   - compeff = compressor efficiency (optional, dimensionless)
            %   - T2 = gas outlet temperature (optional, K)
            %
            %   Heuristics for compressor efficiencies is based on:
            %   Turton et al. (5th Ed.)
            %
            %   Example implementation:
            %   [comppower, compeff, T2] = dsg.sizecompressor(1e5,1,4,323.15,1.5,1.4,0.96)
            %   returns
            %   comppower = 127.4, compeff = 0.8167, T2 = 354.4
            
            if nargin == 6
                Z = 1;
            elseif nargin == 7
                Z = varargin{1};
            end
            
            if P2/P1 > 4
                warning("Compression ratio > 4 is too large - " + ...
                    "check that outlet temperature is not too high! " + ...
                    "Nevertheless continuing calculation...")
            elseif P2/P1 < 1
                error("Outlet pressure smaller than inlet pressure!")
            end
            
            m = m / 3600; % convert kg/h to kg/s
            k = cp/cv;
            a = (k-1)/k;
            power = (m*Z*dsg.R*T1) * ((P2/P1)^a - 1) / a; % useful power
            power = power / 1000; % convert Pa to kPa
            
            compeff = interp1([1,1.5,2,3,6,10], ...
                [0.65-eps,0.65,0.75,0.8,0.85,0.85+eps], P2/P1);
                    
            comppower = power / compeff;
            
            T2 = T1 * (P2/P1)^a;
            
            if T2 > 273.15 + 200
                warning("Gas outlet temperature too high! " + ...
                    "Consider reducing compression ratio P2/P1! " + ...
                    "Nevertheless continuing calculation...")
            end
        end
        
        function [pumppower, pumpeff] = sizepump(Q,dP,varargin)
            % SIZEPUMP Conducts pump sizing by determining required pump power
            %   [pumppower, pumpeff] = SIZEPUMP(Q,dP,rho,pumpeff)
            %   calculates power required to operate a pump based on its
            %   flow rate and pressure differential (discharge - suction
            %   pressure) , where:
            %   - Q = volumetric flow rate through pump (m^3/h)
            %   - dP = pressure differential (kPa)
            %   - rho = stream density (kg/m^3) (optional, default value is
            %   the density of water at room conditions, 1000 kg/m^3)
            %   - pumpeff = pump efficiency calculated (optional, if not
            %   supplied then it will be estimated from heuristics)
            %
            %   Returns:
            %   - pumppower = required pump power (kW)
            %   - pumpeff = pump efficiency (dimensionless) (if pumpeff is
            %   supplied the same value will be returned, otherwise it will
            %   be estimated from heuristics)
            %
            %   Heuristics for pump efficiencies are based on (in order of
            %   priority): (1) GPSA Engineering Databook, (2) North
            %   Carolina Cooperative Extension Service, (3) 75% (assumed
            %   constant) from CAPCOST program.
            %
            %   Example implementation:
            %   [pumppower, pumpeff] = dsg.sizepump(35,500)
            %   returns
            %   pumppower = 8.133 and pumpeff = 0.5977
            
            if nargin == 2
                rho = 1000; % assume 1000 kg/m^3
                pumpeff = -1; % create dummy flag to calculate pumpeff later
            elseif nargin == 3
                rho = varargin{1};
                pumpeff = -1; % create dummy flag to calculate pumpeff later
            elseif nargin == 4
                rho = varargin{1};
                pumpeff = varargin{2};
                if pumpeff < 0 || pumpeff > 1
                    pumpeff = -1; % create dummy flag to calculate pumpeff later
                    warning("Supplied pump efficiency is not between 0 & 1! Recalculating based on heuristics!");
                end
            else
                error("Wrong number of inputs! See documentation!")
            end
            
            power = (Q / 3600) * dP; % useful power in kW
            
            if pumpeff < 0 % catch dummy flag to calculate pumpeff
                H = dP / (rho * dsg.g); % required head in m
                H_ft = H * 3.281; % required head in ft
                Q_gpm = Q * 4.403; % flowrate in gal/min (gpm)
                
                if 50 <= H_ft && H_ft <= 300 && 100 <= Q_gpm && Q_gpm <= 1000
                    pumpeff = ([80, -0.2855, 3.78e-4, -2.38e-7, 5.39e-4, -6.39e-7, 4e-10] * ...
                        [1, H_ft, H_ft*Q_gpm, H_ft*Q_gpm^2, H_ft^2, H_ft^2*Q_gpm, H_ft^2*Q_gpm^2]') / 100;
                else
                    % Maximum useful power for centrifugal pumps = 300 kW
                    pumpeff = interp1([0,2,5,10,30,55,300], ...
                        [0.55-eps,0.55:0.05:0.75,0.75+eps], power);
                end
            end
            
            pumppower = power / pumpeff;
            
        end
        
        function [A, F] = sizeHE_heater(mc,cpc,Tcin,Tcout,Thin,Thout,U,vargin)
            % SIZEHE_HEATER Conducts shell-and-tube heat exchanger sizing, where cold process stream is heated, by determining required heat exchange area
            %   [A, F] = SIZEHE_HEATER(mc,cpc,Tcin,Tcout,cph,Thin,Thout,U,F,Ns)
            %   calculates heat exchange area of shell-and-tube heat
            %   exchanger (counterflow arrangement) for heating a cold
            %   process stream, where:
            %   - mc = cold stream mass flow rate (kg/h)
            %   - cpc = heat capacity of cold stream % J/(kg.K)
            %   - Tcin = cold stream inlet temperature (degC)
            %   - Tcout = cold stream outlet temperature (degC)
            %   - Thin = hot stream inlet temperature (degC)
            %   - Thout = hot stream outlet temperature (degC)
            %   - U = heat transfer coefficient (W/(m^2.degC))
            %   - F = user-specified correction factor (if not specified, F
            %   will be calculated)
            %   - Ns = number of shell passes (default value is 1)
            %   - Ta = temperature approach
            %
            %   Returns:
            %   - A = required heat exchange area
            %   - F = correction factor (optional output - if F is not
            %   specified in input, F will be calculated)
            %
            %   Example implementation:
            %   [A, F] = dsg.sizeHE_heater(31715,3246,89,101,160,156,850)
            %   returns
            %   A = 6.43, F = 0.998
            
            if nargin == 7
                F = -1; % dummy flag to calculate F later
                Ns = 1; % default one shell pass
            elseif nargin == 8
                F = vargin{1};
                Ns = 1; % default one shell pass
            elseif nargin == 9
                F = vargin{1};
                Ns = vargin{2};
            end
            
            % Now, check if temperatures are feasible
            if Thout > Thin
                error("Hot stream outlet cannot be hotter than inlet!")
            elseif Tcout < Tcin
                error("Cold stream outlet cannot be colder than inlet!")
            elseif Tcout > Thout
                warning("Potential temperature cross - Cold stream outlet is hotter than hot stream inlet! Nevertheless continuing with calculations...")
            elseif Thout - Tcin < dsg.Ta
                error("Minimum temperature not fulfilled for hot outlet / cold inlet side!")
            elseif Thin - Tcout < dsg.Ta
                error("Minimum temperature not fulfilled for hot inlet / cold outlet side!")
            end
            
            mc = mc / 3600; % convert kg/h to kg/s
            
            Q = mc * cpc * (Tcout - Tcin); % calculate heat transfer rate
            if (Thout - Tcin) == (Thin - Tcout)
                LMTD = Thout - Tcin;
            else
                LMTD = ((Thin-Tcout) - (Thout-Tcin)) / ...
                    log((Thin-Tcout) / (Thout-Tcin));
            end
            
            if F < 0
                R = (Thin - Thout) / (Tcout - Tcin);
                P = (Tcout - Tcin) / (Thin - Tcin);
                if R == 1
                    W = (Ns-Ns*P) / (Ns-Ns*P+P);
                    F = (sqrt(2)*(1-W)/W) / ...
                        log((W/(1-W) + 1/sqrt(2)) / (W/(1-W) - 1/sqrt(2)));
                else
                    W = ((1-P*R) / (1-P)) ^ (1/Ns);
                    S = sqrt(R^2+1) / (R-1);
                    F = S*log(W) / log((1+W-S+S*W) / (1+W+S-S*W));
                end
            end
            
            A = Q / (U * F * LMTD);
            
        end
        
        function [A, F] = sizeHE_cooler(mh,cph,Thin,Thout,Tcin,Tcout,U,vargin)
            % SIZEHE_COOLER Conducts shell-and-tube heat exchanger sizing, where hot process stream is cooled, by determining required heat exchange area
            %   [A, F] = SIZEHE_COOLER(mh,cph,Thin,Thout,cpc,Tcin,Tcout,U,F,Ns)
            %   calculates heat exchange area of shell-and-tube heat
            %   exchanger (counterflow arrangement) for cooling a hot
            %   process stream, where:
            %   - mh = hot stream mass flow rate (kg/h)
            %   - cph = heat capacity of hot stream % J/(kg.K)
            %   - Thin = hot stream inlet temperature (degC)
            %   - Thout = hot stream outlet temperature (degC)
            %   - Tcin = cold stream inlet temperature (degC)
            %   - Tcout = cold stream outlet temperature (degC)
            %   - U = heat transfer coefficient (W/(m^2.degC))
            %   - F = user-specified correction factor (if not specified, F
            %   will be calculated)
            %   - Ns = number of shell passes (default value is 1)
            %   - Ta = temperature approach
            %
            %   Returns:
            %   - A = required heat exchange area
            %   - F = correction factor (optional output - if F is not
            %   specified in input, F will be calculated)
            %
            %   Example implementation:
            %   [A, F] = dsg.sizeHE_cooler(31715,3246,89,60,5,10,850)
            %   returns
            %   A = 14.8, F = 0.9944
            
            if nargin == 7
                F = -1; % dummy flag to calculate F later
                Ns = 1; % default one shell pass
            elseif nargin == 8
                F = vargin{1};
                Ns = 1; % default one shell pass
            elseif nargin == 9
                F = vargin{1};
                Ns = vargin{2};
            end
            
            % Now, check if temperatures are feasible
            if Thout > Thin
                error("Hot stream outlet cannot be hotter than inlet!")
            elseif Tcout < Tcin
                error("Cold stream outlet cannot be colder than inlet!")
            elseif Tcout > Thout
                warning("Potential temperature cross - Cold stream outlet is hotter than hot stream inlet! Nevertheless continuing with calculations...")
            elseif Thout - Tcin < dsg.Ta
                error("Minimum temperature not fulfilled for hot outlet / cold inlet side!")
            elseif Thin - Tcout < dsg.Ta
                error("Minimum temperature not fulfilled for hot inlet / cold outlet side!")
            end
            
            mh = mh / 3600; % convert kg/h to kg/s
            
            Q = mh * cph * (Thin - Thout); % calculate heat transfer rate
            if (Thout - Tcin) == (Thin - Tcout)
                LMTD = Thout - Tcin;
            else
                LMTD = ((Thin-Tcout) - (Thout-Tcin)) / ...
                    log((Thin-Tcout) / (Thout-Tcin));
            end
            
            if F < 0
                R = (Thin - Thout) / (Tcout - Tcin);
                P = (Tcout - Tcin) / (Thin - Tcin);
                if R == 1
                    W = (Ns-Ns*P) / (Ns-Ns*P+P);
                    F = (sqrt(2)*(1-W)/W) / ...
                        log((W/(1-W) + 1/sqrt(2)) / (W/(1-W) - 1/sqrt(2)));
                else
                    W = ((1-P*R) / (1-P)) ^ (1/Ns);
                    S = sqrt(R^2+1) / (R-1);
                    F = S*log(W) / (log(1+W-S+S*W) / log(1+W+S-S*W));
                end
            end
            
            A = Q / (U * F * LMTD);
            
        end

    end
end