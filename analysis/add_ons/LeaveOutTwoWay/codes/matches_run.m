%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This m-file is used to provide an example on how to compute leave 
%out estimates on a two-way fixed effects model using a test sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %External Packages

%Note: cd is LeaveOutTwoWay                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'~/matlab_bgl/'); %note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
path(path,'codes'); %this contains the main LeaveOut Routines.
warning('off', 'all') 
%cd codes;    
%path(path,'CMG'); %this contains the main LeaveOut Routines.
%MakeCMG;
%cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %PARALLEL ENVIRONMENT
                  
%Note: The user should decide which set-up is most suitable given
%      her own configuration. Make sure to delete the pool once
%      estimation has been carried out.                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
%delete(gcp('nocreate'))
%pool=parpool(14,'IdleTimeout', Inf);
pool=parpool(32,'IdleTimeout', Inf);
%pool=parpool(str2num(getenv('SLURM_CPUS_PER_TASK')),'IdleTimeout', Inf);
%pool = parpool('dcs', 64);
%[a b]=system('sinfo -p low -o "%C"');
%cores=str2num([b(19) b(20)]);
%cores=min(cores,64);
%pool = parpool('dcs', cores);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %CALL
                        
%Note: This is where you should specify where is the input .csv file
%      and how you would like to call the log-file created by the
%      leave out functions. 
%      
%      The user should also specify where to save and name:
%      1. Log File
%      2. Saved Results (will be in .csv)
        

%Make sure that the input .csv file is sorted by worker id and year
%(xtset id year in Stata). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=3;
table=zeros(12,P*1);

for balanced_original_sample = 1 : 1
colonna=1;    
for places = 8 : 8    
for movers_all = 1 : 1     
for type_leave_out = 1 : 1
for type_SE = 1 : 1 

     
if places == 1
namecity='Rovigo'
end
if places == 2
namecity='Belluno'
end
if places == 4
namecity='Venezia'
end
if places == 5
namecity='Padova'
end
if places == 6
namecity='Verona'
end
if places == 7
namecity='Treviso'
end
if places == 8
namecity='Vicenza'
end
if places == 3
namecity='Rovigo_Belluno'
end


%Input File
if balanced_original_sample == 1
bilanciato='BALANCED'    
namesrc=['/scratch/public/leave_out/' namecity '_1996balanced2001_T_greater_2score_test_data.csv']; %where original data is'
end
if balanced_original_sample == 2
bilanciato='UNBALANCED'      
namesrc=['/scratch/public/leave_out/' namecity '_1996unbalanced2001_T_greater_2score_test_data..csv']; %where original data is'
end

%Log
placelog='logs/';
namelog=[namecity '_test_SEs_balanced_panel'];

%Saved File
nameFile=namelog;
placeFile='results/';
filename=[placeFile nameFile];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=importdata(namesrc);
id=data(:,1);
firmid=data(:,2);
year=data(:,3);
y=data(:,4);
clear data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %RUN LEAVE-OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Options (See Description in leave_out_COMPLETE)
if type_leave_out == 1
leave_out_level='obs';
end
if type_leave_out == 2
leave_out_level='matches';
end
if type_leave_out == 3
leave_out_level='workers';
end
andrews_estimates=0;
eigen_diagno=0;
if movers_all == 1 
restrict_movers=1;
moversonly='MOVERS_ONLY'
end
if movers_all == 2 
restrict_movers=0;
moversonly='EVERYONE'
end
resid_controls=0;
controls=[]; %equivalent to no controls
do_SE=0;
subsample_llr_fit=type_SE;
type_of_algorithm='exact';
epsilon=0.005;
eigen_fast=0;
do_montecarlo=0;


%Log File
namelog=['SCORE_TEST' bilanciato '_' moversonly '_' namecity '_leave_out_level_' leave_out_level '_subsample_llr_fit_' num2str(subsample_llr_fit) 'do_standard_errors_' num2str(do_SE) '_algorithm_' type_of_algorithm]
logname=[placelog namelog  '.log'];
system(['rm ' logname])
diary(logname)

%Run
if type_leave_out <= 2
         sigma2_psi = match_test(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE,type_of_algorithm,epsilon,filename)    
end

diary off
end
end
colonna=colonna+1;
end
end
out=table;
s=['results/SE_results_batch_clustering_Hutchison_10000' bilanciato '.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
end

%Close
%delete(pool)