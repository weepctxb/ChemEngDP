# ChemEngDP
MATLAB / Python (U/C) code for mechanical design, sizing &amp; capex/opex calculations

The intention is to use these quickly as "black-box" functions for the purposes for TAC optimization in the second half of the detail design of the main equipment (e.g. PFR length, selection of pumps/compressors etc.)

- ***QUICK START***:

  - Download the required .m files under `\lib` directly to your project path; or

  - Download the entire `\lib` to path, then right click \lib > Add to Path > Selected Folders in MATLAB

To use, download the required .m files into the same folder as your existing code, and call the following functions depending on your situation:

Quick tip: Always use `help` to check the required units to avoid unit conversion errors!

------------------------------------------------

## Mechanical Design - dsg

`help dsg` - retrieves help documentation on functions, inputs and arguments

Functions:

- dsg.designhorzpres - perform entire mechanical design for horizontal pressure vessels

- dsg.designvertpres - perform entire mechanical design for vertical pressure vessels

- dsg.designvac - perform entire mechanical design for vacuum vessels

- dsg.sizecompressor - conducts compressor sizing

- dsg.sizepump - conducts pump sizing

- dsg.sizeHE_heater - conducts heat exchanger sizing for heating stream

- dsg.sizeHE_cooler - conducts heat exchanger sizing for cooling stream

To call intermediate functions (e.g. calculate shell thickness, calculate max. allowable stress, calculate wind allowance etc.), refer to documentation (call `help dsg`).

------------------------------------------------

## CAPEX Calculation - capex

`help capex` - retrieves help documentation on functions, inputs and arguments

Constants:

- capex.CEPCIyy - retrieves annual CEPCI index for 20yy (where yy = 01, 18 or 19)

- capex.USSGyy - retrieves USD:SGD forex rate for 20yy year-average

- capex.CPIzzyy - retrieves country zz's consumer price index (CPI) for 20yy year-average (where zz = 'SG' or 'US')

  - Reference year for CAPCOST = 2001 (as of Turton et al. 5th Ed.)

  - Reference year for utilities cost = 2018 (as of Turton et al. 5th Ed.)

Functions:

- capex.eqptpurcost - calculates purchased equipment cost

- capex.pressurefac - calculates pressure factor

- capex.baremodfac - calculates bare module factor

- capex.baremodcost - calculates bare module cost

- capex.totmodcost - calculates total module cost

- capex.grasscost - calculates grassroots cost

To do: Import relevant data on constants for MOC and pressure factors as needed, only if there are significant benefits of doing so.

------------------------------------------------

## OPEX Calculation - opex

`help opex` - retrieves help documentation on functions, inputs and arguments

Constants:

- opex.CEPCIyy - retrieves annual CEPCI index for 20yy (where yy = 01, 18 or 19)

- opex.USSGyy - retrieves USD:SGD forex rate for 20yy year-average

- opex.CPIzzyy - retrieves country zz's consumer price index (CPI) for 20yy year-average (where zz = 'SG' or 'US')

  - Reference year for CAPCOST = 2001 (as of Turton et al. 5th Ed.)

  - Reference year for utilities cost = 2018 (as of Turton et al. 5th Ed.)

- opex.SF - retrieves stream factor for plant operation

- opex.runtime - retrieves operational runtime per annum

- opex.shiftdur - retrieves duration of one workshift

- opex.shiftperweek - retrieves number of shifts per year

- opex.yearww - retrieves number of work weeks per year after leave entitlements

- opex.utilxxx - utility cost for xxx (USD/GJ basis), where xxx =

  - "HPS" for high-pressure steam (41 barg, 254 degC)

  - "MPS" for medium-pressure steam (10 barg, 184 degC)

  - "LPS" for low-pressure steam (5 barg, 160 degC)

  - "CW" for cooling water (30-45 degC)

  - "ChW" for chilled water (5 degC)

  - "LTR" for low temperature refrigerant (-20 degC)

  - "VLTR" for very low temperature refrigerant (-50 degC)

  - "elec" for electricity (110-440 V)

Functions:

- opex.operatorspershift - calculates # operators per shift

- opex.labourcost - calculates annualised labour cost

- opex.costofutil - calculates cost of utilities

- opex.costofmanfc - calculates all components of annualised total cost of manufacture (COM)