********************************************************************************
Date: 02/04/2022
Email: rsaggio@mail.ubc.ca
********************************************************************************

*DESCRIPTION*     
                    
********************************************************************************
This folder contains all the computer programs needed to replicate the 
results of Di Addario, Kline, Saggio, Sølvsten (DKSS henceforth)
forthcoming in the Journal of Econometrics.

There are two separate parts:

BUILD: This is performed in Stata,  see the master do-file: "build/main.do"  
ANALYSIS: This is performed both in Matlab and Stata, see the master m-file: "analysis/main_INVIND.m" and "analysis/main.do" 
********************************************************************************


*BUILD*     
                    
********************************************************************************
The original raw INPS-INVIND data can be obtained by contacting the Bank of Italy
Research Department.

Three input files should be stored inside the folder build/src:

- "giantNEW_INVIND": contains all the job spells data from the INPS-INVIND data from 1990-2015.
- "matr_bitalia_1990_2013NEW_INVIND": contains information on size, sector and other demographics for the universe of Italian Employers.
- "cerved_smaller": contains information on Value Added for a subset of Italian employers.


The do-file “codes/build_RAKIM.do" combines these files to construct a job-person panel as described in the paper.
The do-file “codes/export_matlab" exports this file into a csv that is going to be read by Matlab.
The do-file “codes/build_AKM" builds additional files needed to compute Table 4 of DKSS.
********************************************************************************

*ANALYSIS (MATLAB)*     
                    
********************************************************************************
Estimation of the DWL model is performed in Matlab. There are two key routines:

1. "kss_rakim_build" finds the equivalent of the "leave-one-out" connected set for the DWL model, see appendix of DKSS.
    If the source data is large (e.g. ~1 million firm effects), then this code can take quite a while to complete (i.e., several days).

2.  "kss_rakim_levels" estimates the DWL model.


Log file obtained after running "analysis/main_INVIND.m" returns all the results
printed in Table 1,2,3,4,5,6,of DKSS as well as Table A3 and Table A4
********************************************************************************

*ANALYSIS (STATA)*     
                    
******************************************************************************** 
Analysis in STATA is performed after having fit the DWL model in Matlab.

The master do-file "analysis/main" computes Figure 1,2,3,4,5,6,7,8, as
well as Table 7 and the Appendix material displayed in Figure A1, A2, A3 and A4 and
Table A1 and Table A2.
********************************************************************************

*EXTERNAL PACKAGES*     
                    
******************************************************************************** 

To run the MATLAB code, the user needs two external packages:

 - matlab_BGL: available at http://www.mathworks.com/matlabcentral/fileexchange/10922
 
 - CMG: Installed via the “MakeCMG.m” command

MATLAB 2018a or higher is required for the codes to run appropriately.

Mac users running a 64-bit machine should install the matlab_BGL package 
available here: https://dgleich.wordpress.com/2010/07/08/matlabbgl-osx-64-bit/
and follow the instructions contained in "installation_mac_KSS" provided in the
replication package for Kline Saggio, Sølvsten (2020, Econometrica):
https://eml.berkeley.edu/~pkline/papers/KSS2020_COMP.zip          
    

     
    