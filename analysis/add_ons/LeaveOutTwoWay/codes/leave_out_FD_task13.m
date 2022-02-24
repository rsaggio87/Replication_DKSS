function [sigma2_psi V_sym V_nosym Q_share C_share_1 C_share_2 C_share_3 C_share_4 C_share_5] = leave_out_FD_task13(y,id,firmid,leave_out_level,leave_2_out,controls,type_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,filename)
%% Author: Raffaele Saggio
%Email: raffaele.saggio@berkeley.edu

%% Version:
% 1.0: Wrote documentation. 06.15.2018.

% 2.0: New version:
%                    - Fixed bug when pcg stagnates when performing JL algorithm.
%                    - Rewrote documentation. 
%                    - Standard errors now account for serial correlation
%                      of the error term within worker.


%% DESCRIPTION
%This function computes the leave out estimates of the variance
%of firm effects in two-way fixed effects model as described in KSS.

%This functions permits computation of the leave-out variance of firm 
%effects by partialling out the worker effects
%and working with the model in first differences (or in stacked First
%Differences when T>2, see below).

%This code also permits computation of the variance of firm effects that
%remains unbiased in the presence of serial correlation in the error term
%within each worker (leave_out_level='workers').

%When working with large dataset it is highly recommended to set the option
%type_algorithm='JL', in order to apply the randomized algorithm routine 
%described in Appendix B of KSS. 
%Also, the user should set the option eigen_fast=1 in order to speed up 
%calculations of the eigenvalues and eigevectors.

%When the input data set has max(T_i)=2, where T_i is the total number of
%person year observations in which we observe a worker, the function works
%with a simple model in First Differences (FD). When max(T_i)>2, we work
%with a model that combines all first differences for a given individual and
%notice that weighted least squares estimates in this Pooled First Differences
%(PFD) model are numerically equivalent to standard within group 
%(Fixed Effects) estimates, provided that each difference for an individual
%is weighted by T_i.

%The code also includes a small montecarlo exercise at the end as described
%in the empirical section of KSS.

%If controls are specifiied by the user the function will work as follows:
%After finding the largest connected set the code estimates a standard AKM
%model with controls. Then it uses the estimated effects on these controls to 
%partial out the effect of these variables on the outcomes. The leave-out
%model will then work with this residualized outcome to speed up
%computation.
%% DESCRIPTION OF THE INPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%y: outcome. Dimensions: N* x 1; N*= # of person-year observations.

%id: worker indicators. Dimensions: N* x 1

%firmid: firm indicators. Dimensions: N* x 1

%leave_out_level: string variable that takes two values:

            %'obs': perform leave-out by leaving a person-year observation
            %out.

            %'workers': perform leave-out by leaving an entire worker's history out.
            
%When leave_out_level='workers', the code automatically adjusts standard 
%errors to account for serial correlation in the error term.
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                    %---NON-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%type_algorithm: 'string'

%This declares which type of algorithm to perform in order to
%compute (Bii,Pii). 

%If type_algorithm='exact', then (Bii,Pii)
%are computed using an exact method which can take quite a bit of time in 
%large datasets. 

%If type_algorithm='JL', then the code uses the Spielman and Srivastava (2008) 
%algorithm which make use of the Johnson-Lindenstrauss lemma. 
%See Appendix B of KSS.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%controls:
%Matrix of controls with dimensions: N* x K. This matrix of controls must
%be appropriately defined by the user ex-ante. For instance, if the user 
%wants to include time effects then the user should include in the matrix 
%'controls' the set of dummy variables associated to a particular year
%effect. If 'controls' is empty, then no controls will be used for
%estimation.


%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                  
%eigen_diagno: Binary. 

% If 1, the code outputs the lindeberg condition and 
% eigenvalue ratio of theorem 1. The code will also output the 
% weak-id confidence intervals using the AM method described in the paper 
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%eigen_fast: Binary. 
%If 1 (and eigen_diagno=1), the code uses a simulation method to estimate 
%the sum of squared eigenvalues of the design matrix.
%It is recommended to set eigen_fast=1 when working with large datasets.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%do_montecarlo: Binary. 
%If 1, the code performs a MC experiment.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%restrict_movers: Binary. 
%If 1, the code performs estimation (and weights the results) using movers
%only.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%filename: string. 
%Used to name saved outputs.
%Default is 'leave_out_FD_estimates';
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-


%% DESCRIPTION OF THE OUTPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%sigma2_psi: leave-out variance of firm effects. 
%V:sampling variance of the leave-out variance of firm effects

%Check Log File for additional results reported by the code. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

no_controls=0;
%% READ
if nargin < 4
error('More arguments needed');
end

if nargin == 4
    no_controls=1;
    controls=ones(size(y,1),1);
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    restrict_movers=0;
    filename='leave_out_FD_estimates';
end

if nargin == 5
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    restrict_movers=0;
    filename='leave_out_FD_estimates';
end

if nargin == 6
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
    restrict_movers=0;
end

if nargin == 7
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
    restrict_movers=0;
end

if nargin == 8
    do_montecarlo=0;
    restrict_movers=0;
    filename='leave_out_FD_estimates';
end

if nargin == 9
    restrict_movers=0;
    filename='leave_out_FD_estimates';
end

if nargin == 10
    filename='leave_out_FD_estimates';
end

if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
end



%% Identify the type of panel (T=2 or T>2).
[~,~,id_norm]=unique(id);
T=accumarray(id_norm,1);
clear id_norm
maxT=max(T);

if maxT==2 % T=2.
    type_estimator='FD';
end

if maxT>2 % General case. Will set model in PFD
    type_estimator='FE';
end


if do_montecarlo==1 && strcmp(type_estimator,'FE')  
    error('Montecarlo can be computed only for T=2 case.');
end



%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp('Listing options')
no_controls
leave_out_level
type_algorithm
type_estimator
eigen_diagno
eigen_fast
do_montecarlo
restrict_movers
filename
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)


%% STEP 1: PRELIMINARIES
%As first step in our analysis, we run estimation of a standard AKM model
%on the original input data. 

%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker

%Find the connected set. Then define dimensions and relevant matrices.
[y,id,firmid,id_old,firmid_old,controls] = connected_set(y,id,firmid,lagfirmid,controls);


%Define
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
L=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
toc
b=pcg(xx,xy,1e-10,1000,L,L');
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
%Do some cleaning of matrices in memory
clear xx xy L xb pe fe ahat ghat F D S Lchol


%%Deresidualized outcome variable (allows the leave one out methodology to run faster)
if no_controls==0
y=y-X(:,N+J:end)*b(N+J:end);
end
clear X b


%% STEP 2: LEAVE ONE OUT CONNECTED SET
%Here we compute the leave out connected set as defined in Appendix B. 
%The input data is represented by the largest connected set. After applying
%the function 'pruning_unbal_v3', the output data will be a connected set
%such that the associated bipartite graph between workers and firms remains
%connected after removing any particular worker from the graph.


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


%Very Important Auxiliaries
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
clear id_mover



s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
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


%% STEP 3: LEAVE TWO OUT CONNECTED SET
tic
	if leave_2_out == 1
	firmid_leave=firmid;
	[id,firmid,y,controls,id_orig,firmid_orig] = pruning_leave2out(id,firmid,y,controls);
	
	%%%Drop stayers with a single person year observation
	T=accumarray(id,1);
	T=T(id);
	sel=T>1;
	y=y(sel,:);
	firmid=firmid(sel,:);
	id=id(sel,:);
	id_orig=id_old(sel,:);
	firmid_orig=firmid_orig(sel,:);
	controls=controls(sel,:);

	%Resetting ids one last time.
	[~,~,n]=unique(firmid);
	firmid=n;
	[~,~,n]=unique(id);
	id=n; 


	%Very Important Auxiliaries
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
	clear id_mover

	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s);
	s=['Info on the leave two out connected set:'];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
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
	end
toc			
%% STEP 4: RESHAPING DATA
%If max(T_i)=2, we simply set the model in FD.
%If max(T_i)>2, we set the model in PFD and use the appropriate weighting.
	
	NT			=	size(y,1);
	F			=	sparse(1:NT,firmid',1);
	J			=	size(F,2);
	count		=	ones(length(id),1);
	gcs 		= 	cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
	


	%Standard FD transformation
	if strcmp(type_estimator,'FD')
		%outcome variable    
		ylag=[NaN; y(1:end-1)];
		ylag(gcs==1)=NaN; %%first obs for each worker
		ydelta=y-ylag;
		ydelta(isnan(ydelta(:,1)),:)=[]; %remove first observation for worker;

		%matrix of assignements
		Flag=[NaN(1,J); F(1:end-1,:)];
		Flag(gcs==1,1)=NaN; %just put it for first firm, we'll remove rows where there is at least one NaN.
		Fdelta=F-Flag;
		sel=~any(isnan(Fdelta),2);
		Fdelta=Fdelta(sel,:);

		%ids
		firmid_delta=[NaN; firmid(1:end-1)];
		firmid_delta(gcs==1)=NaN;
		firmid_delta=firmid_delta(sel);
		firmid_delta_f=firmid;
		firmid_delta_f=firmid_delta_f(sel);
		id_movers=id(sel);
		clear Flag ylag
	end

	%PD transformation
	if strcmp(type_estimator,'FE')
		tic
		[ydelta, Fdelta, Tweight,id_movers,firmid_delta,firmid_delta_f,gcs_delta]= stacked_Fdelta(y,id,firmid,gcs);
		disp('Time to set up the model in Pooled Differences')
		toc
	end

%Design matrix (Laplacian) and preconditioner
L=(Fdelta'*Fdelta); %Laplacian matrix
pfun_ = cmg_sdd(L); %preconditioner for Laplacian matrices.

%% STEP 5: IDENTIFY THE NON ZERO VALUES IN HAT MATRIX
%The next step is probably the most abtstract, yet the most important,
%passage of the code. Leave one out computation requires computation of the
%(block) diagonal of the "hat" matrix, P. This is potentially a
%gigantic matrix. Even if it is block diagonal, some computers might not be
%able to store all its elements. However, thanks to the particular
%structure associated with the Laplacian we know exactly which entries 
%need to be filled in this matrix. 

%If the user perfoms leave-out dropping just a person-year observation,
%then the associated (Bii,Pii) will be a non-zero scalar if the worker
%moved across firms in the corresponding year-year combination. 

%If the user perfoms leave-out dropping the entire history of a worker,
%then the associated (Bii,Pii) will each be a matrix as we have to keep
%track of cross-products across year-year combinations where the worker has
%moved across firms.

tic
if strcmp(type_estimator,'FD')
    index=[(1:size(firmid_delta,1))']; %index of observations in the FD world
    sel=firmid_delta~=firmid_delta_f;
    MOVERS_INDEX=sel;
    elist=[firmid_delta(sel) firmid_delta_f(sel) firmid_delta(sel) firmid_delta_f(sel)]; %keeping track of firm in t, firm in t' for movers. Since here T=2, no need for cross-products when doing leave-out
    rows=index(sel);
    column=index(sel);
    [elist, ~, index_unique]=unique(elist,'rows');
end

if strcmp(type_estimator,'FE') && strcmp(leave_out_level,'obs') 
    index=[(1:size(firmid_delta,1))']; %index of observations in the PFD world
    sel=firmid_delta~=firmid_delta_f;
    elist=[firmid_delta(sel) firmid_delta_f(sel) firmid_delta(sel) firmid_delta_f(sel)]; %keeping track of firm in t, firm in t' for movers. Since we want to leave one obs at the time, no need for cross-products when doing leave-out
    rows=index(sel);
    column=index(sel);
    Tweight=Tweight(sel);
    [elist, ~, index_unique]=unique(elist,'rows');
end

if strcmp(type_estimator,'FE') && strcmp(leave_out_level,'workers')
    [elist,index_unique,rows,column,Tweight]= input_lambda_P_FE_fast(firmid_delta,firmid_delta_f,id_movers,Tweight); %%keeping track of firm in t, firm in t' for movers across periods for a given worker. Now we need also to keep track of cross-products because we leave want to run leave-out leaving the entire worker history of a given individual
end

disp('Time to set up the indexes of the Hat matrix')
toc
disp('Completed Step 4')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%% STEP 5: COMPUTE (Bii,Pii)
%There are two ways to proceed. 

%First way, is to compute exact estimates,
%by parallelizing computation of S_xx^(-1)x_i' across cores.

%Second way uses the SS (2011) algorithm described in Appendix B. Notice
%that everything applies of that Appendix because we still have S_xx has a
%Laplacian matrix in this context.

%We make use of the cmg solver for both steps. The user can fix a seed for 
%replication purposes. 


%Inputs
tol=1e-5; %tol for pcg
epsilon=0.01; %rules the tol for random projection (essentially how many simulations to take).

%Calculate (Bii,Pii)
[Pii, Bii] = eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_);
Pii=Pii(index_unique);
Bii=Bii(index_unique);

%need reweight differences if running FE via PD.
if strcmp(type_estimator,'FE')
    Pii=Pii./Tweight; 
    Bii=Bii./Tweight;
end

%Lambda P
Lambda_P=sparse(rows,column,Pii,size(Fdelta,1),size(Fdelta,1));
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.
%Lambda B
Lambda_B=sparse(rows,column,Bii,size(Fdelta,1),size(Fdelta,1));
Lambda_B=Lambda_B+triu(Lambda_B,1)'; %make it symmetric.
disp('Completed Step 5')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
%% STEP 6: ESTIMATION
%Step 6A: AKM (Plug-in) estimates of the variance of firm effects.
%Step 6B: Verify that AKM and PDF give back same answer.
%Step 6C: Leave out with associated confidence interval (q=0 and q=1).

%% STEP 6A: AKM ESTIMATION on LEAVE ONE OUT LARGEST CONNECTED SET
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix
D=sparse(1:NT,id',1);
X=[D,F*S];
xx=X'*X;
xy=X'*y;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
N=size(D,2);
ghat=b(N+1:N+J-1);
fe=F*S*ghat;
sigma_2_psi_AKM=var(fe);
clear X b xx Lchol ahat ghat S
%% STEP 6B: VERIFY model in PFD gives back same estimates as Fixed Effects.
xy=Fdelta'*ydelta;
b=pcg(L,xy,1e-10,1000,pfun_);
psi_hat_norm=b;
disp('checking - must report 0')
abs(var(fe)-var(F*b)) %you can also check that each single firm effect (when consistently normalized across methods) are numerically equivalent.
r=ydelta-Fdelta*b;
clear xy

%% STEP 6C: LEAVE OUT ESTIMATES
%Leave Out Residuals
Ndelta=size(Fdelta,1);
I_Lambda_P=(speye(Ndelta,Ndelta)-Lambda_P);
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[eta_h, flag]=pcg(I_Lambda_P,r,1e-5,1000,L_P,L_P');

%Auxiliary
A_b=F'*(fe-mean(fe));
[W_to_use, my_first_part] = construc_W_FD(ydelta,Fdelta,L,pfun_,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);

tic
sigma2_psi=leave_out_estimation_two_way_FD(ydelta,Fdelta,F,L,pfun_,sigma_2_psi_AKM,Lambda_B,Lambda_P,I_Lambda_P,L_P,eta_h,W_to_use,my_first_part);


disp('Time to Compute Leave Out of Firm Effects')
toc

%% STEP 7: REPORTING
s=['Results:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(size(F,2))];
disp(s);
s=['# of Person Year Observations: ' num2str(NT)];
disp(s);
s=['Max Leverage' num2str(max(diag(Lambda_P)))];
disp(s);
s=['-*-*-*-*-*-*Variance of Firm Effects-*-*-*-*-*-*'];
disp(s)
s=['AKM: ' num2str(sigma_2_psi_AKM)];
disp(s)
s=['Leave One Out: ' num2str(sigma2_psi)];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
	
	
%Brute force C
tic
C = brute_force_C(Fdelta,F);
s =[filename '_FD_NEW'];
save(s)	
%load(s)
disp('time to brute force C')
toc

%Second piece 
[V_sym V_nosym Q_share C_share_1 C_share_2 C_share_3 C_share_4 C_share_5] = paths_network(id_movers,firmid_delta,firmid_delta_f,ydelta,C);
V_nosym 											 					  = V_nosym/NT^2;
V_sym 												  					  = V_sym/NT^2;


%%%%%%%%%%%%%%%%%%%%%%Montecarlo
	
%fit the sigmas first
sigma_i 	= ydelta.*eta_h; 
SIZE		= diag(F'*F);
SIZE		= log(SIZE);
F1			= F(gcs==1,:);
F2			= F(gcs==2,:);
SIZE1 		= F1*SIZE;
SIZE2		= F2*SIZE;
xdata		= [ones(Ndelta,1) diag(Lambda_P) diag(Lambda_B) SIZE1 SIZE2];
fun 		= @(x,xdata)exp(x(1)*xdata(:,1)+x(2)*xdata(:,2)+x(3)*xdata(:,3)+x(4)*xdata(:,4) + x(5)*xdata(:,5));
x0 			= [0,0,0,0,0];
[b] 		= lsqcurvefit(fun,x0,xdata,sigma_i)
stop
sigma_i		= fun(b,xdata);
figure
histogram(sigma_i(xdata(:,2)>0))
saveas(gcf,[ s 'HISTOGRAM.png'])

%Fill zeros
NSIM			= 1000;
theta			= zeros(NSIM,1);
V_sym_SIM		= zeros(NSIM,1);
V_nosym_SIM		= zeros(NSIM,1);
coverage_SYM	= zeros(NSIM,1);
coverage_NOSYM	= zeros(NSIM,1);

%Create the DGP
psi_true	= psi_hat_norm*sqrt(sigma2_psi/sigma_2_psi_AKM);
sigma2_true = var(F*psi_true);
CC			= C;

for s=1:NSIM
	ydelta						  = Fdelta*psi_true+sqrt(sigma_i).*randn(Ndelta,1);
	
	%Run KSS
	xy                        	  = Fdelta'*ydelta;
	b						  	  = pcg(L,xy,1e-10,1000,pfun_);
	r						      = ydelta-Fdelta*b;
	[eta_h, flag]			      = pcg(I_Lambda_P,r,1e-5,1000,L_P,L_P');
	fe						      = F*b;
	theta(s)				      = var(fe)-(ydelta'*Lambda_B*eta_h)/(NT-1);
	
	%Get the SE
	[V_sym V_nosym] 			  = paths_network(id_movers,firmid_delta,firmid_delta_f,ydelta,CC);
	
	%Normalize
	V_nosym_SIM(s) 				  = sqrt(V_nosym/NT^2);
	V_sym_SIM(s) 				  = sqrt(V_sym/NT^2);
	
	%COVERAGE
	UB							  = theta(s)+1.96*V_sym_SIM(s);
	LB							  = theta(s)-1.96*V_sym_SIM(s);
	coverage_SYM(s)				  = (sigma2_true>=LB && sigma2_true<= UB);
	
	UB							  = theta(s)+1.96*V_nosym_SIM(s);
	LB							  = theta(s)-1.96*V_nosym_SIM(s);
	coverage_NOSYM(s)			  = (sigma2_true>=LB && sigma2_true<= UB);
	
	s
	
end
sel  							  = (imag(V_sym_SIM)>0) + (imag(V_nosym_SIM)>0);
sel								  = sel>0;
s=['Share of draws where I got imaginary SEs'];
disp(s);
mean(sel)


theta			= theta(~sel);
V_sym_SIM		= V_sym_SIM(~sel);
V_nosym_SIM		= V_nosym_SIM(~sel);
coverage_SYM	= coverage_SYM(~sel);
coverage_NOSYM	= coverage_NOSYM(~sel);

oracle_UB=theta+1.96*std(theta);
oracle_LB=theta-1.96*std(theta);
oracle_coverage=(sigma2_true>=oracle_LB).*(sigma2_true<=oracle_UB);
 

s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['True Value of Variance of Firm Effects: ' num2str(sigma2_true)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of Leave one Out Variance Estimator: ' num2str(mean(theta))];
disp(s);
s=['Std of Estimated Variance of Firm Effects: ' num2str(std(theta))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of standard error estimator (Symmetry): '    num2str(mean(V_sym_SIM))];
disp(s);
s=['Expected Value of standard error estimator (No-Symmetry): ' num2str(mean(V_nosym_SIM))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['p1 p10 p25 p50 p75 p90 p99 of mathcal{V}%:']
disp(s)
quantile(V_sym_SIM,[0.01 0.10 0.25 0.50 0.75 0.90 0.99])
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Coverage Rate of oracle estimator 95%: ' num2str(mean(oracle_coverage))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (SYMMETRIC) 95%: ' num2str(mean(coverage_SYM))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (NON-SYMMETRIC) 95%: ' num2str(mean(coverage_NOSYM))];
disp(s);

s=[filename 'MCCC_FD_NEW'];
save(s)


end

