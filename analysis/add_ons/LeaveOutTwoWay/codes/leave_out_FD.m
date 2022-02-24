 function [sigma2_psi, V]  = leave_out_FD(y,id,firmid,leave_out_level,leave_2_out,controls,type_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,epsilon,Q_use,filename)
%% Author: Raffaele Saggio
%Email: raffaele.saggio@berkeley.edu

%% Version:
% 1.0: Wrote documentation. 06.15.2018.

% 2.0:
%                    - Fixed bug when pcg stagnates when performing JL algorithm.
%                    - Rewrote documentation. 
%                    - Standard errors now account for serial correlation
%                      of the error term within worker.

% 3.0: New version after Revision of KSS for ECTA:
%                    - Added option to compute leave 2 out connected set
%                    - New non-parametric inference part based on "split" sample idea. 


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
            
%leave_2_out: binary, if 1 computes estimates on the leave 2 out connected set. We suggest to set this option to 1 if interested in making inference statements of the variance of firm effects                                     
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

%V_symm:sampling variance of the leave-out variance of firm effects imposing that the error term is symmetric

%V_nosymm: sampling variance of the leave-out variance of firm effects without imposing that the error term is symmetric

%Note: Computation of the sampling variance is based on a random selection of an individual travelling along a given edge of the network. Fix seed to replicate results.

%Check Log File for additional results reported by the code. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
no_controls=0;

%% READ
if nargin < 5
error('More arguments needed');
end

if nargin == 5
    no_controls=1;
    controls=ones(size(y,1),1);
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    restrict_movers=0;
    epsilon = 0.005;
    filename='leave_out_FD_estimates';
    Q_use=1;
end

if nargin == 6
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    restrict_movers=0;
     epsilon = 0.005;
    filename='leave_out_FD_estimates';
    Q_use=1;
end

if nargin == 7
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
    epsilon = 0.005;
    restrict_movers=0;
    Q_use=1;
end

if nargin == 8
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
     epsilon = 0.005;
    restrict_movers=0;
    Q_use=1;
end

if nargin == 9
    do_montecarlo=0;
    restrict_movers=0;
        epsilon = 0.005;
    filename='leave_out_FD_estimates';
    Q_use=1;
end

if nargin == 10
    restrict_movers=0;
        epsilon = 0.005;
    filename='leave_out_FD_estimates';
    Q_use=1;
end

if nargin == 11
    epsilon = 0.005;
    filename='leave_out_FD_estimates';
    Q_use=1;
end

if nargin == 12
    filename='leave_out_FD_estimates';
    Q_use=1;
end

if nargin == 13
    filename='leave_out_FD_estimates';
    Q_use=1;
end


if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
end

if size(controls,2)>0
    no_controls=0;
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


if nargout>1 && strcmp(leave_out_level,'workers')  
    warning('Inference is only heteroskedatic robust');
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


%%Deresidualized outcome variable (allows the leave one out methodology to run faster)
if no_controls==0
NT			=   size(firmid,1);
F			=	sparse(1:NT,firmid',1);
D			=	sparse(1:NT,id',1);
J			=  size(F,2);
N			= size(D,2);
S			=	speye(J-1);
S			=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
X			=   [D,F*S,controls];
xx			=   X'*X;
L			=   ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy			=   X'*y;
b			=  pcg(xx,xy,1e-10,1000,L,L');
psii		= b(N+1:J-1);
y			=  y-X(:,N+J:end)*b(N+J:end);
end
clear X xx xy b


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
		id_old_movers=id_old(sel);
		clear Flag ylag
		controls_delta=controls(sel,:);
			
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

%% STEP 4: DIAGNOSTICS FOR ESTIMATING THE VARIANCE OF FIRM EFFECTS
%This part is only computed provided that the option 'eigen_diagno' is
%turned on. In this part, we calculate the squared eigenvalue ratio and
%Lindeberg conditions that constitute the key conditions to verify the
%validity of Theorem 1. Lindebergc condition assumes q=1.

if eigen_diagno==1  

    Faux=Fdelta;
    Faux(:,J)=[];
    xx=Faux'*Faux;
    Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
    
	Sxxdot=[xx sparse(J-1,1); sparse(1,J-1) 1];
	a=F;
	a(:,J)=[];
	abar=sum(a,1)'/(sqrt(NT));
	S_aa=a'*a;
	ainv=pcg(xx,abar,1e-10,1000,Lchol,Lchol');
	Adot=[S_aa -S_aa*ainv;abar' -abar'*ainv];
	[Q, lambda_eig] = eigs(Adot,Sxxdot,3);
	lambda_eig=diag(lambda_eig);
	
	if eigen_fast==0 %To calculate sum of squared eigenvalues (here using exact method)
        Degree=(F'*F);
        Dsqrt=sqrt(Degree);
        Dsqrt_inv=Dsqrt^(-1);
        normL=Dsqrt_inv*L*Dsqrt_inv;
        EIG = eig(full(normL));
        EIG = sort(EIG);
        EIG=EIG(2:end); %Laplacian has always first eigenvalue equal to 0.
        EIG=(1./EIG); %eigenvalues for the inverse.
        EIG=EIG.^2;
        SUM_EIG=sum(EIG);
        clear normL Degree Dsqrt_inv Dsqrt
            for pp=1:3
                EIG_NORM(pp,1)=EIG(pp)/SUM_EIG;
            end
            
    end
    
    if eigen_fast==1 %To calculate sum of squared eigenvalues (via simulations)
        SUM_EIG=trace_Atilde_sqr_FD(Fdelta,F,L,pfun_); 
            for pp=1:3
                EIG_NORM(pp,1)=(lambda_eig(pp,1)^2)/SUM_EIG;
            end    
    end
	
	%Q_use = sum(EIG_NORM>0.10);
	
	for qqq=1:Q_use
		v(:,qqq)=Q(1:end-1,qqq)-ainv*Q(end,qqq);
	end
	
	lambda_1=lambda_eig(1:Q_use);
	x1bar=Faux*v;
    norm=(sum(x1bar.^2)).^(0.5);
    x1bar=x1bar./norm;
    disp('checking, must be identity matrix')
    x1bar'*x1bar
    clear Q xx Faux 
    
    disp('Time to find Eigenvalues and Eigenvectors to compute Lindeberg')
    toc

    
end
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)


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

%% STEP 6A: AKM ESTIMATION
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
clear X b xx Lchol ahat 

%% STEP 6B: VERIFY model in PFD gives back same estimates as Fixed Effects.
xy=Fdelta'*ydelta;
b=pcg(L,xy,1e-10,1000,pfun_);
psi_hat_norm=b;
disp('checking - must report 0')
abs(var(fe)-var(F*b)) %you can also check that each single firm effect (when consistently normalized across methods) are numerically equivalent.
r=ydelta-Fdelta*b;
clear xy
XD=Fdelta;
XD(:,end)=[];
xx=XD'*XD;
xy=XD'*ydelta;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
ydelta_pred=XD*b;
clear XD Lchol xy xx 


%% STEP 6C: LEAVE OUT 
%Leave Out Residuals
Ndelta=size(Fdelta,1);
I_Lambda_P=(speye(Ndelta,Ndelta)-Lambda_P);
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[eta_h, flag]=pcg(I_Lambda_P,r,1e-5,1000,L_P,L_P');

%Auxiliary
A_b=F'*(fe-mean(fe));
[W_to_use, my_first_part] = construc_W_FD(ydelta,Fdelta,L,pfun_,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);

%Done!
sigma2_psi=sigma_2_psi_AKM-(1/NT)*ydelta'*Lambda_B*eta_h;


movers_sel		= firmid_delta~=firmid_delta_f; %focus on movers only
%Standard Error Estimation
if nargout > 1
		
	
		vec				= C_build_FD(ydelta,Fdelta,L,pfun_,Lambda_B,I_Lambda_P,L_P,F);
		%CCC			= brute_force_C(Fdelta,F,lambda_1,x1bar);
		
		
		if eigen_diagno == 1
			Lambda_B2=Lambda_B;
			
			for qqq=1:Q_use
				norm_fact		= lambda_1(qqq)*(x1bar(:,qqq).^2);
				norm_fact		= spdiags(norm_fact,0,size(ydelta,1),size(ydelta,1));
				Lambda_B2		= Lambda_B2-norm_fact;
			end
			
			vec2			= C_build_FD(ydelta,Fdelta,L,pfun_,Lambda_B2,I_Lambda_P,L_P,F);

			for qqq=1:Q_use
			vec2			= vec2 - (lambda_1(qqq)*x1bar(:,qqq))*sum(x1bar(:,qqq).*ydelta); %Update calculation of vec2=B_2*y-0.5(Lambda_B2*eta_h+(I-P)(1-Lambda_P)^(-1)*Lambda_B2)*y = B*y - lambda_1*w1*w1'*y-0.5(Lambda_B2*eta_h+(I-P)(1-Lambda_P)^(-1)*Lambda_B2)*y; where Lambda_B2 = Lambda_B - diag(lambda_1 wi1^2).
			end
			sigma_i			= ydelta.*eta_h;
		
		
		[V, share_bp, COV_R1, A_i2_in_Aj1, A_i1_in_Aj2, A_i2_in_Aj2, Ppath_1, Ppath_2, no_A2, prob_pairs, j_in_Ai1, j_in_Ai2, i_in_Aj1, i_in_Aj2, Cij, D, magical_number, index_movers,index_movers_i,index_movers_2,index_movers_i_2,taxonomy,Dq,magical_number_q] = paths_network_v3(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,x1bar,lambda_1,vec2,sigma_i,SUM_EIG);
		V
		COV_R1
		
		[V share_bp COV_R1] 		  = MC_fixed_FD_v3(ydelta,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,ydelta,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,Cij,D,vec,id_movers(movers_sel),movers_sel,index_movers,index_movers_i,index_movers_2,index_movers_i_2,Ppath_1,Ppath_2,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,no_A2,magical_number,prob_pairs,taxonomy,vec2,Dq,magical_number_q,sigma_i,x1bar,Lambda_B2,lambda_1);
		V
		COV_R1
		
		
		%[V_sym, V_nosym, share_bp, COV_R1_sym, Ppath_1, Ppath_2, no_A2, prob_pairs, j_in_Ai1, j_in_Ai2, i_in_Aj1, i_in_Aj2, Cij, D, magical_number, index_movers,index_movers_i,index_movers_2,index_movers_i_2,index_circle1,index_circle2,focus_circle1_i,focus_circle2_i,focus_circle1_j,focus_circle2_j,PGIGA_CIRCLE1,PGIGA_CIRCLE2,Dq,magical_number_q] = paths_network_v3_old(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,x1bar,lambda_1,vec2,sigma_i,SUM_EIG);
		
		%disp('Relative distance b/w Symmetric Estimator and current estimator')
		%abs((V_sym-V)/V_sym)
		
		%disp('Relative distance b/w non-Symmetric Estimator and current estimator')
		%abs((V_nosym-V)/V_nosym)
		
		%disp('Relative distance b/w non-Symmetric Estimator and current estimator in \SIGMA')
		%abs((COV_R1-COV_R1_sym)./COV_R1_sym)
		
		COV_R1(1:Q_use,Q_use+1) = COV_R1(1:Q_use,Q_use+1)/NT;
		COV_R1(Q_use+1,1:Q_use)	= COV_R1(Q_use+1,1:Q_use)/NT;
		COV_R1(Q_use+1,Q_use+1) = COV_R1(Q_use+1,Q_use+1)/(NT^2);
		b_1						= x1bar'*ydelta;
		theta_1					= sigma2_psi;
		for qqq=1:Q_use
				theta_1			= theta_1-(lambda_1(qqq)/(NT))*(b_1(qqq)^2-COV_R1(qqq,qqq));	
		end
		
		end
		
		if eigen_diagno == 0
		[V, share_bp, COV_R1] = paths_network_v3(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B);
		V
		COV_R1
		end
			
%Time-out for normalization		
		V 			= V/NT^2;
	
end


%% STEP 7: REPORTING
s=['Results: Leave One Out Largest Connected Set'];
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

if eigen_diagno==1
	for pp=1:1
		if pp == 1
		title='Diagnostics on Variance of Firm Effects';
		end    
		if pp == 2
		title='Diagnostics on Variance of Person Effects';
		end
		if pp == 3
		title='Diagnostics on CoVariance of Person, Firm Effects';
		end
		s=['*********************' title '*********************'];
		disp(s);
		s=['ratio of eigenvalues: '];
		disp(s)
		EIG_NORM(1:3,pp)
		s=['Weak Lindeberg Condition: ' num2str(max(x1bar.^2))];
		disp(s)
		s=['Sum of squared eigenvalues: ' num2str(SUM_EIG/NT^2)];
		disp(s)
		if nargout > 1
		s=['Variance of b_1:'];
		disp(s);
		COV_R1(1:Q_use,1:Q_use)
		s=['Variance of theta_1: ' num2str(COV_R1(Q_use+1,Q_use+1))];
		disp(s);
		s=['Covariance of (b_1,theta_1): '];
		disp(s);
		COV_R1(Q_use+1,1:Q_use)
		end
	end
end
s=['******************************************'];
disp(s);

if nargout>1
for pp=1:1
	if pp == 1
	title='Inference on Variance of Firm Effects';
	end    
	if pp == 2
	title='Inference on Variance of Person Effects';
	end
	if pp == 3
	title='Inference on CoVariance of Person, Firm Effects';
	end
	s=['*********************' title '*********************'];
	disp(s);
	s=['SE under q=0: ' num2str(sqrt(V))];
	disp(s);
	s=['Condition for consistency: ' num2str(share_bp)];
	disp(s);;
	UB=sigma2_psi+1.96*sqrt(V);
	LB=sigma2_psi-1.96*sqrt(V);
	s=['CI under q=0: ' num2str(sigma2_psi-1.96*sqrt(V)) '  '  '  ' num2str(sigma2_psi+1.96*sqrt(V))];
	disp(s);;   
	if eigen_diagno==1		
		[UB_AM,LB_AM,Curvature]=AM_CI(NT,lambda_1,COV_R1,b_1,theta_1)
		
		if Q_use  == 1
		gamma_sq=((lambda_1^2/NT^2)*(COV_R1(1,1)^2))/(COV_R1(2,2));
		[UB_AM,LB_AM,Curvature]=AM_CI_old(NT,lambda_1,gamma_sq,COV_R1,b_1,theta_1)
		end
		s=['CI under q=1: ' num2str(LB_AM) '  '  '  ' num2str(UB_AM)];
		disp(s); 
		s=['Curvature: ' num2str(Curvature)];
		disp(s); 
	end
end
end

%%%%%%%%%%%%%%%%%%%%%%Save File
tic
s=[filename 'FD_completed'];
save(s)
toc

if do_montecarlo == 1

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
fun 		= @(x,xdata)exp(xdata(:,1)*x(1)+x(2)*xdata(:,2)+x(3)*xdata(:,3)+x(4)*xdata(:,4)+x(5)*xdata(:,5));
x0 			= [0,0,0,0,0];
b 			= lsqcurvefit(fun,x0,xdata,sigma_i)
sigma_i_true= fun(b,xdata);
psi_true	= psi_hat_norm*sqrt(sigma2_psi/sigma_2_psi_AKM);


%Output the key parameters for the parametric bootstrap ---will only read this if we are working with pooled Veneto sample in leave 1 out.
out			= [id_old_movers,Fdelta*psi_true,sigma_i_true,ydelta];
out			= full(out);
out			= out(movers_sel,:);
s			= [filename 'DGP_values.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
clear out

%Read the parameters of the DGP.
namesrc				= ['results/All_VenetoDGP_CREATION1_2019_7_10_7_35_40DGP_values.csv'];
data				= importdata(namesrc);
id_pool				= data(:,1);
sel					= ismember(id_pool,id_old_movers);
if sum(sel) ~= Nmovers
error('something went wrong when merging the parameters of the DGP from the pooled sample')
end
sigma_i_true		= data(sel,3);
psi_true_diff		= data(sel,2);
data				= [id_pool(sel) id_old_movers(movers_sel) id_movers(movers_sel) sigma_i_true psi_true_diff data(sel,4) ydelta(movers_sel)];
num2str(data(1:100,:))
psi_true_diff		= sparse(id_movers(movers_sel),1,psi_true_diff,Ndelta,1);
sigma_i_true		= sparse(id_movers(movers_sel),1,sigma_i_true,Ndelta,1);
xy					= Fdelta'*psi_true_diff;
psi_hat_norm		= pcg(L,xy,1e-10,1000,pfun_);
psi_true			= psi_hat_norm;
sigma2_true 		= var(F*psi_true)
clear data

out			= [id_old_movers,Fdelta*psi_true,sigma_i_true,ydelta];
out			= full(out);
out			= out(movers_sel,:);
s			= [filename 'DIAGNOSTICS.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
clear out

if 1 == 1
%Fill zeros
NSIM			= 1000;
theta			= zeros(NSIM,1);
theta_homo		= zeros(NSIM,1);
theta_AKM		= zeros(NSIM,1);
V_SIM			= zeros(NSIM,1);
coverage		= zeros(NSIM,1);
coverage_AM		= zeros(NSIM,1);
coverage_oracle	= zeros(NSIM,1);
R2				= zeros(NSIM,1);
Curvature		= zeros(NSIM,1);
COV_R1_sim		= zeros(2,2,NSIM);
bar_Beta_1_simul= zeros(NSIM,1);
theta_2_simul	= zeros(NSIM,1);
lunghezza_q0	= zeros(NSIM,1);
lunghezza_q1	= zeros(NSIM,1);
lunghezza_oracle= zeros(NSIM,1);
correlation_sim = zeros(NSIM,1);
Curvature 		= zeros(NSIM,1);
taxnm			= zeros(NSIM,1);


	
%This is to smooth the parallel process
if 1 == 1
tic
D				= parallel.pool.Constant(D);
Cij				= [];
index_movers    = parallel.pool.Constant(index_movers);
index_movers_i  = parallel.pool.Constant(index_movers_i);
index_movers_2  = parallel.pool.Constant(index_movers_2);
index_movers_i_2= parallel.pool.Constant(index_movers_i_2);
Ppath_1			= parallel.pool.Constant(Ppath_1);
Ppath_2			= parallel.pool.Constant(Ppath_2);
no_A2			= parallel.pool.Constant(no_A2);
prob_pairs    	= parallel.pool.Constant(prob_pairs);
magical_number  = parallel.pool.Constant(magical_number);
j_in_Ai1		= parallel.pool.Constant(j_in_Ai1);
j_in_Ai2		= parallel.pool.Constant(j_in_Ai2);
i_in_Aj1		= parallel.pool.Constant(i_in_Aj1);
i_in_Aj2		= parallel.pool.Constant(i_in_Aj2);
Dq				= parallel.pool.Constant(Dq);			
magical_number_q= parallel.pool.Constant(magical_number_q);
A_i2_in_Aj1		= parallel.pool.Constant(A_i2_in_Aj1);
A_i1_in_Aj2		= parallel.pool.Constant(A_i1_in_Aj2);
A_i2_in_Aj2		= parallel.pool.Constant(A_i2_in_Aj2);
taxonomy		= parallel.pool.Constant(taxonomy);

L				= parallel.pool.Constant(L);
Fdelta			= parallel.pool.Constant(Fdelta);
F				= parallel.pool.Constant(F);
Lambda_B		= parallel.pool.Constant(Lambda_B);
I_Lambda_P		= parallel.pool.Constant(I_Lambda_P);
Lambda_P		= parallel.pool.Constant(Lambda_P);
Lambda_B2		= parallel.pool.Constant(Lambda_B2);

disp('Time to send objects in parallel space')
toc
end
tic

 
parfor s=1:NSIM
	ydelta						  = Fdelta.Value*psi_true+(3/5)*sqrt(sigma_i_true).*trnd(5,Ndelta,1);
	R2(s)						  = var(Fdelta.Value*psi_true)/var(ydelta);
	vec 						  = C_build_FD(ydelta,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P,F.Value);
	vec2						  = C_build_FD(ydelta,Fdelta.Value,L.Value,pfun_,Lambda_B2.Value,I_Lambda_P.Value,L_P,F.Value);
	vec2						  = vec2 - (lambda_1*x1bar)*sum(x1bar.*ydelta); %Update calculation of C2*y=B_2*y-0.5(Lambda_B2*eta_h+(I-P)(1-Lambda_P)^(-1)*Lambda_B2)*y = B*y - lambda_1*w1*w1'*y-0.5(Lambda_B2*eta_h+(I-P)(1-Lambda_P)^(-1)*Lambda_B2)*y; where Lambda_B2 = Lambda_B - diag(lambda_1 wi1^2).
	
	%Run KSS-Andrews
	xy                        	  = Fdelta.Value'*ydelta;
	[b flag]					  = pcg(L.Value,xy,1e-10,1000,pfun_);
	r						      = ydelta-Fdelta.Value*b;
	[eta_h, flag]			      = pcg(I_Lambda_P.Value,r,1e-5,1000,L_P,L_P');
	fe						      = F.Value*b;
	theta(s)				      = var(fe)-(ydelta'*Lambda_B.Value*eta_h)/(NT-1);
	theta_AKM(s)				  = var(fe);
	MSE							  = (sum(r.^2))/(size(ydelta,1)-size(L.Value,2)-1);
	theta_homo(s)				  = var(fe)-(trace(Lambda_B.Value)*MSE)/(NT-1);
	sigma_i						  = ydelta.*eta_h;
	
	%Get the SE
	[V taxnm(s) COV_R1] 		  = MC_fixed_FD_v3(ydelta,A_i2_in_Aj1.Value,A_i1_in_Aj2.Value,A_i2_in_Aj2.Value,ydelta,Fdelta.Value,L.Value,F.Value,tol,epsilon,type_algorithm,pfun_,Lambda_P.Value,Lambda_B.Value,Cij,D.Value,vec,id_movers(movers_sel),movers_sel,index_movers.Value,index_movers_i.Value,index_movers_2.Value,index_movers_i_2.Value,Ppath_1.Value,Ppath_2.Value,j_in_Ai1.Value,j_in_Ai2.Value,i_in_Aj1.Value,i_in_Aj2.Value,no_A2.Value,magical_number.Value,prob_pairs.Value,taxonomy.Value,vec2,Dq.Value,magical_number_q.Value,sigma_i,x1bar,Lambda_B2.Value,lambda_1);
	
	%Normalize
	V_SIM(s) 				  	  = sqrt(V/NT^2);
	COV_R1(1,2) 				  = COV_R1(1,2)/NT;
	COV_R1(2,1) 				  = COV_R1(1,2);
	COV_R1(2,2) 				  = COV_R1(2,2)/(NT^2);
	COV_R1_sim(:,:,s)			  = COV_R1;
		
	%Finish up terms for weak ids
	b_1							  = sum(x1bar.*ydelta);
	gamma_sq					  = ((lambda_1^2/NT^2)*(COV_R1(1,1)^2))/(COV_R1(2,2));		
	Fstatistic					  = (b_1^2)/COV_R1(1,1);
	theta_1						  = theta(s)-(lambda_1/(NT))*(b_1^2-COV_R1(1,1));
	correlation_sim(s)			  = COV_R1(1,2)/(sqrt(COV_R1(2,2))*sqrt(COV_R1(1,1)));
	
	%Coverage
	UB							  = theta(s)+1.96*V_SIM(s);
	LB							  = theta(s)-1.96*V_SIM(s);
	coverage(s)				  	  = (sigma2_true>=LB && sigma2_true<= UB);
	lunghezza_q0(s)				  = UB-LB;
	
	%Final stuff for oracle
	bar_Beta_1_simul(s)			  = b_1;
	theta_2_simul(s)			  = theta_1;
	
	disp('done')
	
	
end
toc
s=[filename 'MCCC_FD_NEW_V3'];
save(s)

%Calculate oracle
COV_R1_oracle					  = cov(bar_Beta_1_simul,theta_2_simul);
correlation_oracle				  = corr(bar_Beta_1_simul,theta_2_simul);
gamma_sq_oracle					  = ((lambda_1^2/NT^2)*(COV_R1_oracle(1,1)^2))/(COV_R1_oracle(2,2));
oracle_coverage_fancy			  = zeros(NSIM,1);
Curvature_true					  = zeros(NSIM,1);
parfor s=1:NSIM
    [UB,LB,Curvature_true(s)] 	  = AM_CI(NT,lambda_1,COV_R1_oracle,bar_Beta_1_simul(s),theta_2_simul(s));
    oracle_coverage_fancy(s) 	  = (sigma2_true>=LB).*(sigma2_true<=UB);
    lunghezza_oracle(s)			  = UB-LB;
    
    [UB,LB,Curvature(s)]	  	  = AM_CI(NT,lambda_1,COV_R1_sim(:,:,s),bar_Beta_1_simul(s),theta_2_simul(s));
	coverage_AM(s)			      = (sigma2_true>=LB && sigma2_true<= UB);
	lunghezza_q1(s)				  = UB-LB;
end


%Imaginary problems
s=['Tabulate cases'];
tabulate(taxnm)

oracle_UB=theta+1.96*std(theta);
oracle_LB=theta-1.96*std(theta);
oracle_coverage=(sigma2_true>=oracle_LB).*(sigma2_true<=oracle_UB);


s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['R2 of the model: ' num2str(mean(R2))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['True Value of Variance of Firm Effects: ' num2str(sigma2_true)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of KSS: ' num2str(mean(theta))];
disp(s);
s=['Std of KSS: ' num2str(std(theta))];
disp(s);
s=['Expected Value of Andrews: ' num2str(mean(theta_homo))];
disp(s);
s=['Std of Andrews: ' num2str(std(theta_homo))];
disp(s);
s=['Expected Value of AKM: ' num2str(mean(theta_AKM))];
disp(s);
s=['Std of AKM: ' num2str(std(theta_AKM))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of standard error estimator: '    num2str(mean(V_SIM))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Std Deviation of standard error estimator: '    num2str(std(V_SIM))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Lenght CI under q=0: '    num2str(mean(lunghezza_q0))];
disp(s);
s=['Lenght CI under q=1: ' num2str(mean(lunghezza_q1))];
disp(s);
s=['Lenght CI under q=1 - ORACLE: ' num2str(mean(lunghezza_oracle))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['p1 p10 p25 p50 p75 p90 p99 of mathcal{V}%:']
disp(s)
quantile(V_SIM,[0.01 0.10 0.25 0.50 0.75 0.90 0.99])
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Coverage Rate of oracle estimator 95%: ' num2str(mean(oracle_coverage))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (AM-oracle) 95%: ' num2str(mean(oracle_coverage_fancy))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference 95%: ' num2str(mean(coverage))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (AM) 95%: ' num2str(mean(coverage_AM))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['True Curvature: ' num2str(mean(Curvature_true))];
disp(s);
s=['Average Curvature: ' num2str(mean(Curvature))];
disp(s);
s=['Std Curvature: ' num2str(std(Curvature))];
disp(s);
s=['p1 p10 p25 p50 p75 p90 p99 of Curvature%:']
disp(s)
quantile(Curvature,[0.01 0.10 0.25 0.50 0.75 0.90 0.99])

s=['Oracle SIGMA:'];
disp(s);
num2str(COV_R1_oracle)

s=['Mean Estimated SIGMA:'];
disp(s);
num2str(mean(COV_R1_sim,3))

s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
end
end
end
