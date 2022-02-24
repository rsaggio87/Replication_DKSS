function [slope] = OD_FD_shrink(y,id,firmid,controls,type_algorithm,simulated,filename)

%This function should old be used with a dataset where T=2.


%This objective of this function is to evaluate the performance of a shrunk
%firm effects. It does that by performing an out-of-sample exercise where
%we try to see how the wage changes of workers belonging to a validation
%sample can be predicted on the basis of changes in firm effects estimated
%with a training sample.

%The shrunk firm effects are computed at the firm level. 


if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
end

if size(controls,2)>0
    no_controls=0;
    resid_controls=1;
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
    error('We can only run this when T = 2 currently')
end


%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp('Listing options')
no_controls
type_algorithm
type_estimator
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

%Leave One Out Largest Connected Set
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave one out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic
[y,firmid,id,~,~,controls,~] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);
disp('Time to find leave one out largest connected set')
toc
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%Auxiliary 
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
firm_size =accumarray(firmid,1);
firm_size(end)=[];

%We will be focusing on movers only
if 1 == 1 
    y=y(movers,:);
    firmid=firmid(movers,:);
    id=id(movers,:);
    controls=controls(movers,:);
end

%reset ids
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n;

%Training - Orig Sample indicator
if 0 == 1
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
firmid_delta=[NaN; firmid(1:end-1)];
firmid_delta(gcs==1)=NaN;
firmid_delta=firmid_delta(gcs~=1);
firmid_delta_f=firmid;
firmid_delta_f=firmid_delta_f(gcs~=1);
[~,~,OD_ID]=unique([firmid_delta firmid_delta_f],'rows');
id_movers=id(gcs~=1);
[~,IX]=sort(OD_ID);
OD_ID=OD_ID(IX);
id_movers=id_movers(IX);
count=ones(length(id_movers),1);
gcs  = cell2mat(accumarray(OD_ID,count,[],@(x){cumsum(x)}));
tabulate(gcs)
gcs  = mod(gcs,2);
id_training   = id_movers(gcs == 1);
id_validation = id_movers(gcs == 0);
size(id_training,1)
size(id_validation,1)
size(id_movers,1)
yold=y;
firmid_old=firmid;
id_old=id;

%Training sample
sel = ismember(id,id_training);
y=yold(sel,:);
id=id(sel,:);
firmid=firmid(sel,:);
controls=controls(sel,:);
[~,~,n]=unique(id);
id=n;

%Validation sample
sel = ismember(id,id_validation);
y_valid=yold(sel);
firmid_valid=firmid_old(sel);
id_valid=id_old(sel);
out=[y_valid,id_valid,firmid_valid];
s=[filename 'validation.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
end

s=['Info on the Leave Out Training Sample Sample:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(max(id))];
disp(s);
s=['# of Firms: ' num2str(max(firmid))];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%Partial Out the effects of controls in this data: 
if no_controls == 0
NT=size(y,1);
D=sparse(1:NT,id',1);
N=size(D,2);
F=sparse(1:NT,firmid',1);
J=size(F,2);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
X=[D,F*S,controls];
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=X'*y;
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
y=y-X(:,N+J:end)*b(N+J:end);
end


%% STEP 3: RESHAPING DATA
count=ones(length(id),1);
gcs = cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));


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
    clear Flag ylag
end

%Design matrix (Laplacian) and preconditioner
L=(Fdelta'*Fdelta); %Laplacian matrix
pfun_ = cmg_sdd(L); %preconditioner for Laplacian matrices.


%% STEP 4:  Pii of the AKM model in First Differences; Bii for for variance of firm effects (unweighted)
Ndelta=size(ydelta,1);
index=[(1:size(firmid_delta,1))']; %index of observations in the FD world
sel=firmid_delta~=firmid_delta_f;
elist=[firmid_delta(sel) firmid_delta_f(sel) firmid_delta(sel) firmid_delta_f(sel)]; %keeping track of firm in t, firm in t' for movers. Since here T=2, no need for cross-products when doing leave-out
rows=index(sel);
column=index(sel);
[elist, ~, index_unique]=unique(elist,'rows');

%Calculate (Pii,Bii)
[Pii, Bii] =  eff_res_FD_shrink(elist,Fdelta,L,pfun_);
Pii=Pii(index_unique);
Bii=Bii(index_unique);
%Lambda P
Lambda_P=sparse(rows,column,Pii,size(Fdelta,1),size(Fdelta,1));
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.
%Lambda B
Lambda_B=sparse(rows,column,Bii,size(Fdelta,1),size(Fdelta,1));
Lambda_B=Lambda_B+triu(Lambda_B,1)'; %make it symmetric.

%% STEP 5:  ESTIMATE AKM MODEL IN FIRST DIFFERENCES
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
Fdelta=Fdelta*S;
J=J-1;
L=Fdelta'*Fdelta;
Lchol = ichol(L,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=Fdelta'*ydelta;
[psi,flag]=pcg(L,xy,1e-10,1000,Lchol,Lchol');
psi=psi-mean(psi);
eta=ydelta-Fdelta*psi;
I_Lambda_P=(speye(Ndelta,Ndelta)-Lambda_P);
eta_h=I_Lambda_P\eta; %Leave one out residual
msg = lastwarn ; 
if contains(msg, 'singular') 
            	s=['******************************************'];
                disp(s);
                disp(s); 
				disp('Warning: OLS coefficient not always identified when leaving a particular set of observation out as specified by "leave_out_level"')
				disp('One example where this occurs is when the user asks to run leave out on matches without restricting the analysis to movers only.')
				s=['******************************************'];
                disp(s);
                disp(s); 
				
end
sigma_i = ydelta.*eta_h;
MSE=sum(eta.^2)/(Ndelta-J-1);
%% STEP 5.1:  Estimate of variance
sigma2_KSS  =  var(psi)-(sum(Bii.*sigma_i))/(J-1)
sigma2_homo =  var(psi)-(sum(Bii)*MSE)/(J-1)

%% STEP 6:  Shrinkage
Linv=L^(-1);
[V,~]=eigs(L,2);
Xnew=Fdelta*V;
hBest=1/Ndelta^(1/3);
f = fit(Xnew,sigma_i,'lowess','Normalize','on','Span',hBest);
sigma_predict = feval(f,Xnew);
selnan=isnan(sigma_predict);
sigma_predict(selnan)=sigma_i(selnan);
sel=find(sigma_predict<0);
sigma_predict(sel)=abs(sigma_predict(sel));

%Multivariate Shrinkage - HOMOSKEDASTIC
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,MSE,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_homo.*speye(J)+V);
numeratore       = sigma2_homo.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde        = PI_hat*psi;

%Multivariate Shrinkage - KSS 
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,sigma_i,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_KSS.*speye(J)+V);
numeratore       = sigma2_KSS.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde_kss    = PI_hat*psi;

%Multivariate Shrinkage - KSS Smoothed
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,sigma_predict,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_KSS.*speye(J)+V);
numeratore       = sigma2_KSS.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde_kss_s  = PI_hat*psi;

%Multivariate Shrinkage - KSS Mixed
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,MSE,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_KSS.*speye(J)+V);
numeratore       = sigma2_KSS.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde_kss_M  = PI_hat*psi;


out=[ydelta,Fdelta*psi,Fdelta*psi_tilde,Fdelta*psi_tilde_kss,Fdelta*psi_tilde_kss_s,Fdelta*psi_tilde_kss_M,firmid_delta,firmid_delta_f];
s=[filename 'SHRUNK_training.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);

out=[log(firm_size),psi,psi_tilde,psi_tilde_kss,psi_tilde_kss_s,psi_tilde_kss_M];
s=[filename 'SHRUNK_firm_size.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
%% STEP 7:  Simulation
if simulated == 1
SIMUL = 1000;  

%Truth
sigma2_psi_KSS_true=sigma2_KSS;
[V,~]=eigs(L,2);
Xnew=Fdelta*V;
hBest=1/Ndelta^(1/3);
f = fit(Xnew,sigma_i,'lowess','Normalize','on','Span',hBest);
sigma_predict = feval(f,Xnew);
selnan=isnan(sigma_predict);
sigma_predict(selnan)=sigma_i(selnan);
sel=find(sigma_predict<0);
sigma_predict(sel)=abs(sigma_predict(sel));
%sigma_predict = MSE.*(Fdelta*psi>0)+1.2*MSE.*(Fdelta*psi<=0);
sigma_i_orig=sigma_predict;

%Fill zeros
slope_uncorr=zeros(SIMUL,1);
slope_oracle=zeros(SIMUL,1);
slope_homo=zeros(SIMUL,1);
slope_kss=zeros(SIMUL,1);
slope_kss_s=zeros(SIMUL,1);
slope_kss_M=zeros(SIMUL,1);
estim2=zeros(SIMUL,1);


%Oracle Pi
sigma_i_sparse=sparse((1:Ndelta)',1:Ndelta,sigma_i_orig,Ndelta,Ndelta);
tic
V=Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
toc
denominatore = (sigma2_psi_KSS_true.*speye(J)+V);
numeratore   = sigma2_psi_KSS_true.*speye(J);
PI           = denominatore\numeratore;


slope_size_oracle=zeros(SIMUL,1);
slope_size_uncorr=zeros(SIMUL,1);
slope_size_homo=zeros(SIMUL,1);
slope_size_kss=zeros(SIMUL,1);
slope_size_kss_s=zeros(SIMUL,1);
slope_size_kss_M=zeros(SIMUL,1);

%% SIMULATIONS
parfor s = 1 : SIMUL
    
%Draw from DGP    
psi_true = sqrt(sigma2_psi_KSS_true)*randn(J,1);
psi_true = psi_true- mean(psi_true);
ydelta   = Fdelta*psi_true + sqrt(sigma_i_orig).*randn(size(Fdelta,1),1);
firm_size = 0.4*psi_true   + sqrt(2*sigma2_psi_KSS_true).*randn(size(Fdelta,2),1);

%Estimate based on draw
xy=Fdelta'*ydelta;
[psi,flag]=pcg(L,xy,1e-10,1000,pfun_);
psi=psi-mean(psi);
eta=ydelta-Fdelta*psi;
I_Lambda_P=(speye(Ndelta,Ndelta)-Lambda_P);
eta_h=I_Lambda_P\eta; %Leave one out residual
sigma_i = ydelta.*eta_h;
MSE=sum(eta.^2)/(Ndelta-J-1);
COV=cov(psi,psi_true);
slope_uncorr(s)=COV(1,2)/COV(1,1);

%RUN KSS Correction
sigma2_KSS  =  var(psi)-(sum(Bii.*sigma_i))/(J-1);
sigma2_homo = var(psi)-(sum(Bii)*MSE)/(J-1);
estim(s)=sigma2_KSS;
estim2(s)=sigma2_homo;


%Smoothing the sigmas
hBest=10/Ndelta^(1/3);
f = fit(Xnew,sigma_i,'lowess','Normalize','on','Span',hBest);
sigma_predict = feval(f,Xnew);
selnan=isnan(sigma_predict);
sigma_predict(selnan)=sigma_i(selnan);

%Multivariate Shrinkage - HOMOSKEDASTIC
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,MSE,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_homo.*speye(J)+V);
numeratore       = sigma2_homo.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde        = PI_hat*psi;
COV              = cov(psi_tilde,psi_true);
slope_homo(s)    = COV(1,2)/COV(1,1);

%Multivariate Shrinkage - KSS 
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,sigma_i,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_KSS.*speye(J)+V);
numeratore       = sigma2_KSS.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde_kss    = PI_hat*psi;
COV              = cov(psi_tilde_kss,psi_true);
slope_kss(s)     = COV(1,2)/COV(1,1);

%Multivariate Shrinkage - KSS Smoothed
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,sigma_predict,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_KSS.*speye(J)+V);
numeratore       = sigma2_KSS.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde_kss_s  = PI_hat*psi;
COV              = cov(psi_tilde_kss_s,psi_true);
slope_kss_s(s)   = COV(1,2)/COV(1,1);

%Multivariate Shrinkage - KSS Mixed
sigma_i_sparse   = sparse((1:Ndelta)',1:Ndelta,MSE,Ndelta,Ndelta);
V                = Linv*(Fdelta'*sigma_i_sparse*Fdelta)*Linv;
denominatore     = (sigma2_KSS.*speye(J)+V);
numeratore       = sigma2_KSS.*speye(J);
PI_hat           = denominatore\numeratore;
psi_tilde_kss_M  = PI_hat*psi;
COV              = cov(psi_tilde_kss_M,psi_true);
slope_kss_M(s)   = COV(1,2)/COV(1,1);


%Oracle Shrinkage
psi_oracle=PI*psi;
COV=cov(psi_oracle,psi_true);
slope_oracle(s)=COV(1,2)/COV(1,1);

%Regress firm size
COV=cov(psi_oracle,firm_size);
slope_size_oracle(s)=COV(1,2)/COV(1,1);

COV=cov(psi,firm_size);
slope_size_uncorr(s)=COV(1,2)/COV(1,1);

COV=cov(psi_tilde,firm_size);
slope_size_homo(s)=COV(1,2)/COV(1,1);

COV=cov(psi_tilde_kss,firm_size);
slope_size_kss(s)=COV(1,2)/COV(1,1);

COV=cov(psi_tilde_kss_s,firm_size);
slope_size_kss_s(s)=COV(1,2)/COV(1,1);

COV=cov(psi_tilde_kss_M,firm_size);
slope_size_kss_M(s)=COV(1,2)/COV(1,1);

%ydelta   = Fdelta*psi_true + sqrt(sigma_i_orig).*randn(size(Fdelta,1),1);
end

disp('Average KSS estimate of variance of firm effects (unweighted) - First is KSS, Second is Andrews')
mean(estim)
mean(estim2)
sigma2_psi_KSS_true

slope_uncorr=mean(slope_uncorr)
slope_homo=mean(slope_homo)
slope_kss_s=mean(slope_kss_s)
slope_kss=mean(slope_kss)
slope_kss_M=mean(slope_kss_M)
slope_oracle=mean(slope_oracle)



slope_size_oracle=mean(slope_size_oracle)
slope_size_uncorr=mean(slope_size_uncorr)
slope_size_homo=mean(slope_size_homo)
slope_size_kss=mean(slope_size_kss)
slope_size_kss_s=mean(slope_size_kss_s)
slope_size_kss_M=mean(slope_size_kss_M)

% out=[psi_true,psi,psi_tilde,psi_tilde_kss,psi_tilde_kss_s,psi_tilde_kss_M,psi_oracle];
% s=[filename 'SHRUNK_.csv'];
% dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
% 
% out=[ydelta,Fdelta*psi,Fdelta*psi_tilde,Fdelta*psi_tilde_kss,Fdelta*psi_tilde_kss_s,Fdelta*psi_tilde_kss_M,Fdelta*psi_oracle];
% s=[filename 'SHRUNK_delta.csv'];
% dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
end
end

