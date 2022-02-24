%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
%This m-file provides five examples that illustrate how to conduct inference on 
%linear combinations of regression coefficients. Examples 1-2 conduct inference 
%that is robust to heteroscedasticity and many regressors. Examples 3-5 additionally
%allow errors to be correlated within worker-firm match. 

%All examples below refer to a case where the user in a prior step has 
%performed KSS leave out estimation of a two-way fixed effects model. This 
%is because "lincom_KSS" is essentially a post-estimation command, optimized 
%to run following the command "leave_out_COMPLETE". 
%Nevertheless, "lincom_KSS" can be applied to ANY linear regression model.

%The user should have installed all the files contained in the GitHub
%repository at https://github.com/rsaggio87/LeaveOutTwoWay to properly run
%the code below.

%See the documentation inside "lincom_KSS" for further details. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Testing whether Firm Fixed Effects are different across two regions.

%In this example we will start by loading a person-year file that contains
%wage information from two time periods on two Italian regions. 
%The object is to test whether the firm effects from these two regions are 
%statistically different.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%Paths and other Auxiliaries
path(path,'~/matlab_bgl/'); %note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
path(path,'codes'); %this contains the main LeaveOut Routines.
cd codes;    
path(path,'CMG'); %preoconditioner for Laplacian matrices
warning('off','all')
MakeCMG;
cd ..

%Parpooling
%pool=parpool(24,'IdleTimeout', Inf);
[a b]=system('sinfo -p low -o "%C"');
cores=str2num([b(19) b(20)]);
cores=min(cores,61);
pool = parpool('dcs', cores);

%Log File
placelog='logs/';
namelog='MC_homoskedatic_DGP_one_pred_redraw';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%Name the file to be saved
nameFile=namelog;
placeFile='results/';
filename=[placeFile nameFile];

%Import data 
namesrc='/scratch/public/leave_out/Verona_1996or2001_unbalanced.csv'; %where original data is
data=importdata(namesrc);
id=data(:,1);
firmid=data(:,2);
y=data(:,4);
extra_vars=data(:,5:end); %Column 1: gender (1 for female). %Column 2: Workers' Age

%Run Leave Out Estimation 
leave_out_level='matches';
andrews_estimates=0;
eigen_diagno=0;
restrict_movers=1;
resid_controls=2;
controls=extra_vars; %equivalent to no controls
do_SE=0;
subsample_llr_fit=0;
type_of_algorithm='exact';
epsilon=0.01;


%Run
%[sigma2_psi] = leave_out_COMPLETE(y,id,firmid,leave_out_level,extra_vars,2,0,0,0,restrict_movers,0,'exact',epsilon,filename);

%Load the .csv file created by the function  leave_out_COMPLETE
results=importdata([filename '.csv']);
y=results(:,1);
firmid=results(:,2);
id=results(:,3);
extra_vars=results(:,7:end);
[~,~,match_id]=unique([id firmid],'rows');
clear results

%Load the matrix Lambda_P saved by leave_out_COMPLETE (Note: now this matrix is going to be block-diagonal instead of diagonal as in Example 2) 
Lambda_P=load([filename '_Lambda_P']);
Lambda_P_old=Lambda_P.Lambda_P;
Lambda_P=diag(diag(Lambda_P_old)); %just diagonal matrix for this DGP.
max(max(Lambda_P))

%Now create the X associated with the two-way model (dropping as usual last firm).
NT=size(y,1);
D=sparse(1:NT,id',1); 
F=sparse(1:NT,firmid',1);
N=size(D,2);
J=size(F,2);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];
F=F*S;
X=[D,F]; %Design Matrix.

Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
Z=[extra_vars(:,1)];
labels={'Female Worker'};

%In the code v is always defined as v=(Z'*Z)^(-1)*Z'*Transform and the
%function "lincom_KSS" always adds a constant to the user specified matrix
%"Z". Also recall that beta=(X'*X)^(-1)X'y;

%Therefore in this case v'*beta in Remark 9 of KSS gives back the 4 (gender,
%worker age, log firm size, firm age) t-statistics obtained when projecting
%the firm effects onto a constant, gender dummy, Age, Log Firm Size and Firm Age.

%t_stat_example_3=lincom_KSS(y,X,Z,Transform,match_id,Lambda_P,labels); %cluster the standard errors by match identifiers now.

[coverage_MC] = lincom_KSS_MC(y,X,Z,Transform,[],Lambda_P,labels); %homoskedatic DGP
diary off


%Log File
placelog='logs/';
namelog='MC_cluster_DGP_one_pred_redraw';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

Lambda_P = Lambda_P_old;
[coverage_MC] = lincom_KSS_MC(y,X,Z,Transform,match_id,Lambda_P,labels); %Cluster DGP

%Close
diary off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
if 0 == 1

%Log File
placelog='logs/';
namelog='MC_test_example1';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%How to name the results
nameFile='test_example1';
placeFile='results/';
filename=[placeFile nameFile];

%Import data 
namesrc='src/file_example1.csv'; %where original data is
data=importdata(namesrc);
id=data(:,1); 
firmid=data(:,2);
y=data(:,5);
region=data(:,4); %Region indicator. Value -1 for region 1, Value 1 for region 2.

%Run Leave Out Estimation
leave_out_level='obs';
resid_controls=2; 
andrews_estimates=0;
eigen_diagno=0;
subsample_llr_fit=0;
restrict_movers=0;
do_SE=0;
type_of_algorithm='exact';

sigma2_psi = leave_out_COMPLETE(y,id,firmid,leave_out_level,region,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,[],filename);

%Load data associated with Leave One Out Connected Set. 
results=importdata([filename '.csv']);
y=results(:,1);
firmid=results(:,2);
id=results(:,3);
region=results(:,7:end);
clear results

%Load the matrix of statistical leverages saved by leave_out_COMPLETE
Lambda_P=load([filename '_Lambda_P']);
Lambda_P=Lambda_P.Lambda_P;

%Now create the X associated with the two-way model (dropping as usual last firm).
NT=size(y,1);
D=sparse(1:NT,id',1); 
F=sparse(1:NT,firmid',1);
N=size(D,2);
J=size(F,2);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];
F=F*S;
X=[D,F]; %Design Matrix.

%Test linear combination v'beta where beta=(X'*X)X'*y.
Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
region_dummy=region;
region_dummy(region_dummy==-1)=0;
Z=region_dummy; 
labels={'Region 2 Dummy'};

[coverage_MC] = lincom_KSS_MC(y,X,Z,Transform,[],Lambda_P,labels); %homoskedatic DGP
diary off

end
%Close
diary off