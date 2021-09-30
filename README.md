# SingleMoleculeTransportAssay
modified matlab scripts that originally appear in SPARTAN software (Juette et.al., 2013, Nature Methods) for turnover time calculation of single molecule transport trajectories. 
/[Turnover Time]
make sure all responding trajectories and their corresponding response times are in the same folder
and make sure you are calling 'Turnover Time' function from within that folder
create a folder called 'fits'
select the responding traces, press cancel in the second window
then select the excel sheet file (*.xlsx) with the corresponding response files
OPTIONAL: if you do not want to see the fits you can change create plot parameter in the beginning to 'false'
the function will generate 3 files called
coefficients: this file has all the data for almost all trajectories whether they pass fit criteria or not
rate_constants: this file has 13 parameters written in the header of trajectories that pass the rsquare criteria and many other criteria embedded in the script, take a look at them if you like looking at someoneelse's code, it'll be a small torture but you can do it
rate_constant_not: this file has 13 parameters written in the header of trajectories that DO NOT pass the rsquare criteria
