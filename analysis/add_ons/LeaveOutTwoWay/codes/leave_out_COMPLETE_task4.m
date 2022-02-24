function [R2,adjR2,dof] = leave_out_COMPLETE_task4(QQQ,y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename)
%% Author: Raffaele Saggio
%Email: rsaggio@princeton.edu

%% DESCRIPTION
%This function computes leave out estimates in two-way fixed effects 
%model and conducts inference as described in KSS.

%The mandatory input is person-year dataset that has to be sorted
%by workers' identifiers (id) and year (xtset id year in Stata). The function
%automatically performs computation of the largest connected set and leave
%out connected set. 

%This function can be applied to any two-way fixed effects model
%(student-teachers, patient-doctors, etc) and supports unbalanced data as
%well as T>2 panels.

%We use AKM jargon (workers, firms) when describing the code for
%simplicity.
%% Version:
% 1.0: Wrote documentation. 06.15.2018.

% 1.1:  Speed up computation of eigenvalues/vectors - by avoiding storage of
%       large matrices. 06.18.2018.

% 1.11: Eliminated the option Ndiagno. 06.19.2018.

% 1.12: Read Nargout and compute only parameters asked by user and added
%       option for whether user wants computation of standard error.
%       06.22.2018

% 1.2:  Added options to run fast Local Linear Regression, helpful in large
        %datasets. 01.07.2018 

% 1.25: Fix some bugs arising when setting leave_out_level to either
%       "matches" or "workers"  10.07.2018 

% 1.3:  Dropped stayers that have only one person year observations (for
%       which Pii=1 when estimating the model in levels). 25.07.2018

% 1.32: Improved readability of saved results at the end of the code. 
%       Changed locations of saved results 
%       Added export to .csv results 31.07.2018

% 1.5:  Significantly improved computation of Leave out matrices when there
%       are no controls in the model or the controls have been residualized
%       in a prior step. 
%          
%               We introduced the following changes:
%
%                1. Added CMG routine to speed computation of linear system
%                involving the Laplacian matrix as design matrix. Users need
%                to make sure that the CMG package is installed and placed
%                in the main directory as shown in the GitHub repository.
%
%                2. Read movers-stayers structure to fasten computation of (Bii,Pii).
%
%                In terms of speed, for the test dataset used in "example.m": 
%
%                1. With version 1.32 the code takes 260 seconds to compute (Bii,Pii).
%                2. With version 1.5 the code takes  23 seconds to compute (Bii,Pii).
% 
% 
% 1.51: Better management of large sparse matrices when invoking parfor to 
%       compute using the option parallel.pool.Constant.
%
% 1.52: Added more outputs to the function to simplify possible post-estimation 
%       commands.
%
% 1.55: Added more options to run non-parametric fit.
%
% 2.0:  Added option to approximate (Bii,Pii) using Random Projections methods
%       that build on the Johnson Lindestrauss Lemma - See Appendix B of KSS.
%
%       This especially helpful in massive datasets where exact computation 
%       of (Bii,Pii), even after the improvements introduced from version 1.5, 
%       is close to be prohibitive in terms of computation time.
%          
%       In particular, with this new release we added the following inputs:
%
%                1. type_of_algorithm: This takes two values: "exact" or "JLL".
%
%                    "exact": corresponds to exact computation of (Bii,Pii)
%                     as in version 1.5.          
%                    
%                    "JLL": applies random projection methods to
%                    approximate (Bii,Pii) as detailed in Appendix B of
%                    KSS. 
%                   
%                    Default is "exact".
%                
%                 2. epsilon: this governs the tradeoff b/w speed and unbiasdness 
%                    when estimating (Bii,Pii). Smaller values of epsilon implies 
%                    less bias but more computation time.
%
% 2.05: Made additional improvements for computation of (Bii,Pii) in "eff_res" by looking
%		at unique matches.

% 2.07: Residualize step now computed on leave out largest connected set.
%
% 2.15: Fixed a typo in computation of the standard errors. Improved warning
%       system of the code in case  1-Pii is approximately 0.
%
% 2.17: Residualization step now computed on the leave out connected set.
%
% 2.2:  Code can now handle leave out on observations corresponding to a
%       given match (unique pair of a worker and a firm), i.e. when
%       specifying "leave_out_level=matches"
%
%       - When leaving out a given match, person
%         effects are not identified for stayers. The code will issue a
%         warning if the user is setting "leave_out_level=matches" without
%         setting "restrict_movers=1" and will not output results concerning
%         the variance of person effects.

%       - Note that we added non-mandatory input "year" which reports the 
%         calendar year of a given worker-year observation. 
%         This information is used to conduct the smoothing step 
%         for computation of the standard error, see
%         function "smooth_fn" for details. We highly suggest to include a
%         set of year effects as additional controls for the AKM equation
%         and to set resid_controls=1.

%      -  Computation of the SEs requires knowledge of the trace of a large
%         matrix that is unfeasible to store in memory. We therefore
%         estimate this trace via simulations, similarly as when 
%         "leave_out_level=obs", this time using the Hutchsion trace
%         estimator. See https://arxiv.org/pdf/1308.2475.pdf for an
%         introduction on the Hutchsion trace estimator and the derivation
%         of the bounds used inside the function "hutchison_trace_SE".
%
%      -  We have also augmented the code to take into account the rare
%         instances when the KSS formula for the standard error in the high
%         rank case returns an imaginary root.
%
%      -  See this document for an introduction on how to conduct inference
%         of variance components in the presence of serial correlation
%         within matches: XXXX

% 2.21:   New option: do_SE=2 which reports a conservative standard error of
%         the variance components, robust to both serial correlation and
%         heteroskedasticity. The advantage is that computation of this 
%         standard errors is several orders of magnitudes faster compared 
%         to the case where do_SE=1;
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
%% 					DESCRIPTION OF THE INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%y: outcome. Dimensions: N* x 1; N*= # of person-year observations.
%--
%id: worker indicators. Dimensions: N* x 1
%--
%firmid: firm indicators. Dimensions: N* x 1
%--
%leave_out_level: string variable that takes two values:

    %'obs': perform leave-out by leaving a person-year observation out (default)

    %'matches': perform leave-out by leaving an entire person-firm match out.
   
%When leave_out_level='matches', the standard errors for each variance
%components are robust to serial correlation of the error term within a
%specific match.
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                    %---NON-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%year:
%This is a vector of dimension N* x 1 that should report the calendar year
%of a given person year observation. If this information is provided and
%leave_out_level='matches', then the code will read it to specify a 
%flexible smoothing step by allowing each \hat{\sigma}_i^2 to depend 
%arbitrarly on calendar year as well as allowing dependence within each
%match. This information is not read when leave_out_level='obs', see
%function "smooth_fn" for further details.

%
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%controls:
%Matrix of controls with dimensions: N* x K. This matrix of controls must
%be appropriately defined by the user ex-ante. For instance, if the user 
%wants to include time effects then the user should include in the matrix 
%'controls' the set of dummy variables associated to a particular year
%effect, making sure to avoid potential collinearity issues. 

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%resid_controls: 0 1 2. 

%If 0, the code includes the vector of controls specified by the user
%(provided that is not empty) in estimation of the two-way model.

%If 1 and input "controls" is not empty, then the code partials out 
%the effects of these controls before performing leave out estimates.
%In particular, in Step 1, after performing AKM on the largest connected
%set, we partial out the effect of controls using the corresponding AKM
%estimates.

%If 2 and input "controls" is not empty, the code acts as if the controls
%are not used in estimation of the model but saves the corresponding vector
%of controls in the leave out connected set.

%Default is 0. 

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
%andrews_estimates: Binary.

% If 1, the code reports homoskedastic corrections as described of 
% Andrews (2008). Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
%eigen_diagno: Binary. 

% If 1, the code outputs the lindeberg condition and 
% eigenvalue ratio of theorem 1. The code will also output the 
% weak-id confidence intervals using the AM method described in the paper
% (assuming q=1) provided that do_SE=1.

% Default is 0.  

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%subsample_llr_fit: Discrete. 

%This governs how to compute the local linear regressions (LLR) fit needed
%to derive the standard errors as described in KSS when q=0. 
%See Remark 10 and appendix B of KSS.

%If 0, the code computes LLR (i.e. Lowess model) using the entire data 
%(i.e. both movers and stayers). 

%If 1, and there are no controls in the model, the code distinguishes
%between movers and stayers. For stayers, the assigned fitted value 
%is the average of \hat{\sigma}_i within cells defined in
%terms of unique values of T_i, where T_i is the number of person-year
%observations in which we observe a given worker. For movers, the code runs
%a standard LLR model. 

%If 1, and there are controls in the model, the code runs a stratified
%LLR separately for movers and stayers.

%If 2, and there are no controls in the model, the code creates "Kgrid"
%equally sized bins of Bii and Pii for movers only. 
%The code then fits a weighted LLR model across these "KGrid" x "KGrid" cells 
%weighting by cell size. For stayers, the assigned fitted value 
%is the average of \hat{\sigma}_i within cells defined in
%terms of unique values of T_i, where T_i is the number of person-year
%observations in which we observe a given worker. 

%If 2, and there are controls in the model, the code creates "Kgrid"
%equally sized bins of Bii and Pii for both movers and stayers. 
%We then fit a stratified binned LLR model separately 
%for movers and stayers across these "KGrid" x "KGrid" cells weighting 
%the estimates by cell size. 

%If 3, the code creates "Kgrid" equally sized bins of Bii and Pii
%The code then fits a weighted LLR model across these "KGrid" x "KGrid" 
%cells weighting by cell size. 

%If 4, the code fits a multivariate Nadaraya Watson style of non parametric
%estimator by averaging the value of \hat{\sigma}_i within "KGrid" x "KGrid" 
%cells of (Bii,Pii).

%See the function "llr_fit" for further details. Default value of
%"Kgrid=1000".

%All the technique aboves can also account for the case of serial
%correlation within a match. See code below for further details.  

%Default value is 2. 

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%restrict_movers: Binary. 
%If 1, the code performs leave-out estimation just by focusing on sample
%of movers.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%do_SE: do_SE. 
%If 1, the code provides calculations of the standard errors assuming
%strong identification (q=1). 

%If 2, the code reports a conservarive error assuming again strong
%identification (q=1).

% Default is 1.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%type_of_algorithm: 
%This takes two values: "exact" or "JLL".                    

%   "exact": performs exact computation of (Bii,Pii). 

%   "JLL": perform random projection methods to approximate (Bii,Pii) as 
%   detailed in Appendix B of KSS. 

%In large datasets (say where #of Workers + #of Firms > 50K) we suggest
%setting type_of_algorithm="JLL". As described in Table B1 in large
%datasets JLL gives back idenentical estimates compared to exact methods
%while taking 20-30 times less computation time. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%epsilon: This should be b/w 0 and 1.

%epsilon governs the tradeoff b/w speed and accuracy when estimating 
%(Bii,Pii) using type_of_algorithm=JLL.

%Smaller values of epsilon implies less accuracy but higher computation time.

%Default is 0.1. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%filename: string. 
%Where saved results should be stored and named. Use name like 
%"leave_out_results" and not "leave_out_results.csv"

%Default is 'leave_one_out_estimates';
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

% DESCRIPTION OF THE OUTPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%sigma2_psi: leave-out variance of firm effects.
%SE_sigma2_psi: Standard Error (SE) of leave-out variance of firm effects
%(q=0, high rank case)

%sigma_psi_alpha: leave-out covariance of firm, person effects.
%SE_sigma_psi_alpha: SE of leave-out covariance of firm, person effects.

%sigma2_alpha: leave-out variance of person effects. 
%SE_sigma2_alpha: SE of leave-out variance of person effects.

%y: outcome variable in the leave out connected set.

%id: re-normalized worker identifiers in the leave out connected set.

%firmid: re-normalized firm identifiers in the leave out connected set.

%controls: vector of extra controls in the leave out connected set.

%Pii: (diagonal) statistical leverage associated with two way model.

%Check Log File for additional results reported by the code. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

% SAVED FILES

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%The code saves the following files:

%One .mat file after step 3  with all matrices in memory.
%One .mat file after completing the code with all matrices in memory.
%One .mat file containing the relevant leave-out leverages.
%A   .csv file with location and name specified by the user.

%The .csv saves the following variables belonging 
%to the leave out connected set.

%y: outcome variable.
%firmid: firm identifier. (normalized)
%id: worker identifier. (normalized)
%firmid_old: firm identifier. (as in original input data)
%firmid: firm identifier. (as in original input data)
%controls: vector of controls
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%% READ
no_controls=0;
no_year_provided=0;
if nargin < 4
    error('More arguments needed');
end

if nargin == 4
    no_controls=1;
    controls=ones(size(y,1),1);
    resid_controls=0;
    andrews_estimates=0;
    eigen_diagno=0;
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.1;   
    filename='leave_one_out_estimates';
    year = cell2mat(accumarray(id,ones(size(id,1)),[],@(x){cumsum(x)})); %if year is empty filled it up with consecutive count of person-year observations within each worker employment history.
    no_year_provided=1;
end

if nargin == 5   
    no_controls=1;
    resid_controls=1; 
    controls=ones(size(y,1),1);
    andrews_estimates=0;
    eigen_diagno=0;
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;   
    filename='leave_one_out_estimates';
end

if nargin == 6
    andrews_estimates=0;    
    eigen_diagno=0;
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 7    
    eigen_diagno=0;    
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 8      
    subsample_llr_fit=0;    
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 9     
    restrict_movers=0; 
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 10  
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 11 
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 12 
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 13 
    filename='leave_one_out_estimates';
end


if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
    resid_controls=0;
end

if size(year,2)==0
    no_year_provided=1;
    year = cell2mat(accumarray(id,ones(size(id,1)),[],@(x){cumsum(x)}));
end

if size(controls,2)>0 && resid_controls==2
    no_controls=1;
    resid_controls=0;
end


if resid_controls==1 && no_controls== 1 
    error('cannot residualize if there are no controls specified')    
end

if resid_controls==0 && no_controls== 0 && strcmp(type_of_algorithm,'JLL')
    error('cannot run JLL algorithm on non-Laplacian design matrix. Please set resid_controls=1')    
end

if do_SE==2 && eigen_diagno== 1
    error('No conservative inference in the weak identification case. Please set do_SE=1')    
end


%Read number of outputs
if  nargout==1
    n_of_parameters=1;
end

if nargout==2
    n_of_parameters=2;
end

if nargout>=3
    n_of_parameters=3;
end


%Warnings
if ~strcmp(leave_out_level,'obs') && restrict_movers == 0 
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
disp('WARNING!!!!!!!!!!!!!!!')
disp('Person Effects not identified for stayers (individuals that do not change employers) when leaving more than a single observation out.')
disp('The code will omit all output concerning the variance of person effects.')
disp('Consider setting "restrict_movers=1" if user is interested in the variance of person effects.')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
n_of_parameters=2;
sigma2_alpha=NaN;
SE_sigma2_alpha=NaN;
end

if ~strcmp(leave_out_level,'obs') && do_SE == 1 && no_year_provided==1
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
disp('WARNING!!!!!!!!!!!!!!!')
disp('It is highly recommended to specify a variable that contains the calander year of each person-year observations ')
disp('if the user wants inference robust to serial correlation within a match and potential non-stationarities in the error term.')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
end

if ~strcmp(leave_out_level,'obs') && do_SE == 1 && eigen_diagno==1
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
disp('WARNING!!!!!!!!!!!!!!!')
disp('Currently the code is unable to calculate the confidence intervals using the A-M method described in KSS in the presence of serial correlation error. ')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
end


%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp('Listing options')
n_of_parameters
leave_out_level
no_controls
resid_controls
andrews_estimates
eigen_diagno
subsample_llr_fit
restrict_movers
do_SE
type_of_algorithm
epsilon
filename
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)





if ~strcmp(leave_out_level,'obs') && no_controls == 0 && resid_controls == 0
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
disp('WARNING!!!!!!!!!!!!!!!')
disp('To speed-up computation Leave out on matches will be perfomed in 2 steps. ')
disp('In the first step we residualize the controls')
disp('In the second step we calculate (Bii,Pii) for a model that controls for worker and firm characteristics')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
resid_controls=1;
end

%% STEP 1: PRELIMINARIES
%As first step in our analysis, we run estimation of a standard AKM model
%on the original input data. 

%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker

%Find the connected set. Then define dimensions and relevant matrices.
controls=[controls year];

[y,id,firmid,id_old,firmid_old,controls] = connected_set(y,id,firmid,lagfirmid,controls);

%Define
year=controls(:,end);
controls=controls(:,1:end-1);
NT=size(y,1);
J=max(firmid);
N=max(id);
D=sparse(1:NT,id',1);
F=sparse(1:NT,firmid',1);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
if no_controls==1
X=[D,F*S];
end
if no_controls==0
X=[D,F*S,controls];
end

%Run AKM
disp('Running AKM...')
tic
xx=X'*X;
xy=X'*y;
tic
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
toc
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
toc
ahat=b(1:N);
ghat=b(N+1:N+J-1);
pe=D*ahat;
fe=F*S*ghat;
xb=X*b;
r=y-xb;
dof=NT-size(X,2)-1;
TSS=sum((y-mean(y)).^2);
R2=1-sum(r.^2)/TSS;
adjR2=1-sum(r.^2)/TSS*(NT-1)/dof;


%Some auxiliaries for summary statistics.
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker
stayer=(firmid==lagfirmid);
stayer(gcs==1)=1;
stayer=accumarray(id,stayer);
T=accumarray(id,1);
stayer=T==stayer;
movers=stayer~=1;
movers=movers(id);
id_movers=id(movers);
[~,~,n]=unique(id_movers);
Nmovers=max(n);

%Run AKM
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['STEP 1: AKM Estimates on Largest Connected Set '];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(max(Nmovers))];
disp(s);
s=['# of Firms: ' num2str(size(F,2))];
disp(s);
s=['# of Person Year Observations: ' num2str(sum(diag((F'*F))))];
disp(s);
s=['-*-*-*-*-*-*AKM RESULTS-*-*-*-*-*-*'];
disp(s)
COV=cov(fe,pe);
s=['Variance of Firm Effects: ' num2str(COV(1,1))];
disp(s);
s=['Covariance of Firm and Person Effects: ' num2str(COV(1,2))];
disp(s);
s=['Variance of Person Effects: ' num2str((COV(2,2)))];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str((corr(fe,pe)))];
disp(s);
s=['R2: ' num2str(R2)];
disp(s);
s=['Adj.R2: ' num2str(adjR2)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%Do some cleaning of matrices in memory
clear xx xy L xb pe fe ahat ghat F D S Lchol X b r

%% STEP 2: LEAVE ONE OUT CONNECTED SET
%Here we compute the leave out connected set as defined in Appendix B. 
%The input data is represented by the largest connected set. After applying
%the function 'pruning_unbal_v3', the output data will be a connected set
%such that the associated bipartite graph between workers and firms remains
%connected after removing any particular worker from the graph.
controls=[controls year];

%Focus estimation on movers only?
if restrict_movers==1
sel=movers;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);

%reset ids
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n;
end

%little trick to carry information on year in the leave out connected set
%Leave One Out Largest Connected Set
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave one out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic

[y,firmid,id,id_old,firmid_old,controls] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);
disp('Time to find leave one out largest connected set')
toc
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%Bring back year and the original set of controls
year = controls(:,end);
controls = controls(:,1:end-1);

%%%Drop stayers with a single person year observation
T=accumarray(id,1);
T=T(id);
sel=T>1;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);

%Resetting ids one last time.
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n; 

%Important Auxiliaries
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker
stayer=(firmid==lagfirmid);
stayer(gcs==1)=1;
stayer=accumarray(id,stayer);
T=accumarray(id,1);
stayer=T==stayer;
movers=stayer~=1;
movers=movers(id);
T=T(id);
id_movers=id(movers);
[~,~,n]=unique(id_movers);
Nmovers=max(n);
NT=size(y,1);
D=sparse(1:NT,id',1);
N=size(D,2);
F=sparse(1:NT,firmid',1);
J=size(F,2);
clear id_mover
%Summarize
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(max(firmid))];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)


%Residualize
if resid_controls == 1
   S=speye(J-1);
   S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
   X=[D,F*S,controls];
   xx=X'*X;
   Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
   xy=X'*y;
   b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
   y=y-X(:,N+J:end)*b(N+J:end);
   no_controls=1; 
end

if no_controls==0
   K=size(controls,2);
end

if no_controls==1
    K=0;
end

%% STEP 5: ESTIMATION
%We now conduct estimation on leave-out connected set. 

%These are the steps:

%Step 5.1: Here we perform standard AKM (i.e. plug-in) estimates of the
%variance decomposition parameters.


%Step 5.2: Here we perform leave-out estimates and inference of the
%variance decomposition parameters.

%Step 5.3: Here we perform homoskedasticity corrected estimates of the
%variance decomposition parameters.



%% STEP 5.1: AKM (PLUG-IN)
if K == 0
    S=speye(J-1);
    S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix     
    X=[D,F*S];
    xx=X'*X;
    Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
end

xy=X'*y;
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
ahat=b(1:N);
ghat=b(N+1:N+J-1);
pe=D*ahat;
fe=F*S*ghat;
xb=X*b;
eta=y-xb;
dof=NT-size(X,2)-1;
TSS=sum((y-mean(y)).^2);
R2=1-sum(eta.^2)/TSS
adjR2=1-sum(eta.^2)/TSS*(NT-1)/dof

%% Ridge
X=[ones(NT,1) X];
xx=X'*X;
xx_aux=xx+QQQ*speye(size(xx,2));
Lchol=ichol(xx_aux,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=X'*y;
b = pcg(xx_aux,xy,1e-10,1000,Lchol,Lchol');
xb=X*b;
R2=corr(y,xb)^2

%calculate dof now
aux=xx_aux\X';
NT
K=size(xx,2)
dof=trace(X*aux)

end

