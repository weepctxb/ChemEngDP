% Example 1: Horizontal pressure vessel
% Po = 470 psig, To = 850 degF, Di = 78 in, L = 480 in
mechdesign1 = dsg.designhorzpres(470,850,78,480)

% Example 2: Vertical pressure vessel
% Po = 95.5 psig, To = 150 degF, Di = 120 in, L = 2544 in
mechdesign2 = dsg.designvertpres(95.5,150,120,2544)

% Example 3: Vacuum vessel
% Po = 7.977 psig, To = 257 degF, Di = 168 in, L = 1080 in
mechdesign3 = dsg.designvac(7.977,257,168,1080)

% Example 4: Compressor sizing for required power and outlet temperature
% m = 1e5 kg/h, P1 = 100 kPa, P2 = 300 kPa, T1 = 323.15 K,
% cp = 1.02, cv = 0.72, Z = 0.99
[comppower, compeff, T2] = dsg.sizecompressor(1e5,100,300,323.15,1.02,0.72,0.99)

% Example 5: Pump sizing for required power
% Q = 35 m^3/h, dP = 500 kPa, rho = 1000 kg/m^3 (default)
[pumppower, pumpeff] = dsg.sizepump(35,500)

% Example 6: Heat exchanger sizing for required heat exchange area (for
% heating process stream)
% mc = 31715 kg/h, cpc = 3246 J/(kg.K), Tcin = 89 degC, Tcout = 101 degC,
% Thin = 160 degC, Thout = 156 degC, U = 850 W/(m^2.degC), Ns = 1 (default)
[A1, F1] = dsg.sizeHE_heater(31715,3246,89,101,160,156,850)

% Example 7: Heat exchanger sizing for required heat exchange area (for
% cooling process stream)
% mc = 31715 kg/h, cph = 3246 J/(kg.K), Thin = 89 degC, Thout = 60 degC,
% Tcin = 5 degC, Tcout = 10 degC, U = 850 W/(m^2.degC), Ns = 1 (default)
[A2, F2] = dsg.sizeHE_cooler(31715,3246,89,60,5,10,850)