function [Nmovers, J, share, R2, adjR2, KSS_R2 ,KSS_R2_double_movers,KSS_test] = OD_FD(y,id,firmid,leave_out_level,controls,type_algorithm,eigen_diagno,eigen_fast,NSIM_SE,double_movers,mix_adj,filename)

%% EXTERNAL CALLS
%MakeCMG; %preconditioner for Laplacians matrix: http://www.cs.cmu.edu/~jkoutis/cmg.html
path(path,'~/matlab_bgl/'); %path to the matlabBGL files. note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
do_montecarlo=0;
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
end

if nargin == 5
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
end

if nargin == 6
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
end

if nargin == 7
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
end

if nargin == 8
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
end

if nargin == 9
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

%Leave One Out Largest Connected Set
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave one out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic
[y,firmid,id,id_old,firmid_old,controls,~] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);
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

%movers only
y=y(movers,:);
firmid=firmid(movers,:);
id=id(movers,:);
id_old=id_old(movers,:);
firmid_old=firmid_old(movers,:);
controls=controls(movers,:);

%reset ids
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n;
max(n)

%Select only double moves
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    firmid_delta=[NaN; firmid(1:end-1)];
    firmid_delta(gcs==1)=NaN;
    firmid_delta=firmid_delta(gcs~=1);
    firmid_delta_f=firmid;
    firmid_delta_f=firmid_delta_f(gcs~=1);
    [~,~,matc]=unique([firmid_delta firmid_delta_f],'rows');
    cC=accumarray(matc,1);
    cC=cC(matc);
    sel=cC>=2;
    id_movers=id(gcs~=1);
    id_movers=id_movers(sel);
    sel=ismember(id,id_movers);
    share=mean(sel);
    sel_double=sel;

%OD sample
if double_movers==1
    y=y(sel,:);
    firmid=firmid(sel,:);
    id=id(sel,:);
    controls=controls(sel,:);
    [~,~,n]=unique(firmid);
    firmid=n;
    [~,~,n]=unique(id);
    id=n;
    max(n)
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    [y,id,firmid,~,~,controls] = connected_set(y,id,firmid,lagfirmid,controls);
end

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

s=['Info on the OD largest connected set:'];
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
%% STEP 3: RESHAPING DATA
%If max(T_i)=2, we simply set the model in FD.
%If max(T_i)>2, we set the model in PFD and use the appropriate weighting.

NT=size(y,1);
F=sparse(1:NT,firmid',1);
J=size(F,2);
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
    id_movers=id(sel);
    sel_double=sel_double(sel);
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
%% STEP 4:  Pii of the AKM model 
Ndelta=size(ydelta,1);
index=[(1:size(firmid_delta,1))']; %index of observations in the FD world
sel=firmid_delta~=firmid_delta_f;
elist=[firmid_delta(sel) firmid_delta_f(sel) firmid_delta(sel) firmid_delta_f(sel)]; %keeping track of firm in t, firm in t' for movers. Since here T=2, no need for cross-products when doing leave-out
rows=index(sel);
column=index(sel);
[elist, ~, index_unique]=unique(elist,'rows');
size(elist)
%Inputs
tol=1e-5; %tol for pcg
epsilon=0.01; %rules the tol for random projection (essentially how many simulations to take).
[Pii, ~] = eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_);
Pii=Pii(index_unique);
%Lambda P
Lambda_P=sparse(rows,column,Pii,size(Fdelta,1),size(Fdelta,1));
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.



%% STEP 5: OD Model and associated Bii
[~,~,match_OD]=unique([firmid_delta,firmid_delta_f],'rows');
M=sparse(1:Ndelta,match_OD',1);
disp('Total number of OD effects')
max(match_OD) 
mm=M'*M;
aux=mm\M';
tic
Bii_P=zeros(Ndelta,1);
Bii_1=zeros(Ndelta,1);
Bii=zeros(Ndelta,1);
parfor i=1:Ndelta
    mi=M*aux(:,i);
    xy=Fdelta'*mi;
    [b,flag]=pcg(L,xy,1e-10,1000,pfun_); %Project OD effects on firm effects
    res=mi-Fdelta*b;
    Bii(i)=res'*res;
    Bii_P(i)=var(Fdelta*b)*(Ndelta-1);
    Bii_1(i)=var(mi)*(Ndelta-1);
end

toc

%% STEP 6: AKM
xy=Fdelta'*ydelta;
[psi,flag]=pcg(L,xy,1e-10,1000,pfun_);
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
%% STEP 6: Construct the numerator of the test, report auxiliaries
mu=M\ydelta;
mu=M*mu;
xy=Fdelta'*mu;
[psi,flag]=pcg(L,xy,1e-10,1000,pfun_); %Project OD effects on firm effects
res=mu-Fdelta*psi;

%Mixed Adjustments as suggested by Mikkel.
if mix_adj == 1
    res_od=ydelta-mu;
    T_od=accumarray(match_OD,1);
    T_od=T_od(match_OD);
    Pii_od=1./T_od;
    res_od=res_od./(1-Pii_od);
    sigma_od=ydelta.*res_od;
    sigma_i(T_od>=2)=sigma_od(T_od>=2);    
end
test_num=(sum(res.^2)-sum(Bii.*sigma_i))/Ndelta;




dof=Ndelta-size(Fdelta,2)-1;
TSS=sum((mu-mean(mu)).^2);
R2=1-sum(res.^2)/TSS;
adjR2=1-sum(res.^2)/TSS*(Ndelta-1)/dof;
KSS_R2_num=var(Fdelta*psi)-(1/Ndelta)*sum(Bii_P.*sigma_i);
KSS_R2_den=var(mu)-(1/Ndelta)*sum(Bii_1.*sigma_i);
KSS_R2=KSS_R2_num/KSS_R2_den;
KSS_R2_num=var(Fdelta(sel_double,:)*psi)-(1/Ndelta)*sum(Bii_P(sel_double).*sigma_i(sel_double));
KSS_R2_den=var(mu(sel_double))-(1/Ndelta)*sum(Bii_1(sel_double).*sigma_i(sel_double));
KSS_R2_double_movers=KSS_R2_num/KSS_R2_den
disp('R2 when regressing OD effects on worker, firm dummies')
R2
disp('AdjR2 when regressing OD effects on worker, firm dummies')
adjR2
disp('KSS R2 under the null')
KSS_R2
disp('KSS R2 OD cells with 2 movers or more')
KSS_R2_double_movers


%% STEP 7: Standard errors
if nargout>=8
    if mix_adj == 0
        hBest=1/Ndelta^(1/3);
        f = fit([diag(Lambda_P) Bii],sigma_i,'lowess','Normalize','on','Span',hBest);
        sigma_predict = feval(f,[Pii, Bii]);
        selnan=isnan(sigma_predict);
        sigma_predict(selnan)=sigma_i(selnan);
    end
    if mix_adj == 1
        sel_counts=find(T_od>1);
        %More than 1
        Ncels=sum(T_od>1);
        hBest=1/Ncels^(1/3);
        f = fit([Pii_od(sel_counts) Bii(sel_counts)],sigma_i(sel_counts),'lowess','Normalize','on','Span',hBest);
        sigma_aux = feval(f,[Pii_od(sel_counts) Bii(sel_counts)]);
        selnan=isnan(sigma_aux);
        sigma_aux(selnan)=sigma_i(selnan);
        sigma_predict=sparse(sel_counts,1,sigma_aux,Ndelta,1);
        %1
        sel_counts=find(T_od==1);
        Ncels=sum(T_od>1);
        hBest=1/Ncels^(1/3);
        Pii_HO=diag(Lambda_P);
        f = fit([Pii_HO(sel_counts) Bii(sel_counts)],sigma_i(sel_counts),'lowess','Normalize','on','Span',hBest);
        sigma_aux = feval(f,[Pii_HO(sel_counts) Bii(sel_counts)]);
        selnan=isnan(sigma_aux);
        sigma_aux(selnan)=sigma_i(selnan);
        sigma_predict=sigma_predict+sparse(sel_counts,1,sigma_aux,Ndelta,1);
    end
    aux_SIM=zeros(NSIM_SE,1);
     parfor s=1:NSIM_SE
               v=randn(Ndelta,1).*sqrt(sigma_predict);
               coeff=M\v;
               coeff=M*coeff;
               xy=Fdelta'*coeff;
               [b,flag]=pcg(L,xy,1e-10,1000,pfun_);
               resD=coeff-Fdelta*b;
               subtract=resD'*resD;
               aux=v-coeff;
               aux=aux./(1-Pii);
               aux=v.*aux;
               aux_SIM(s)=subtract-sum(Bii.*aux);
     end
    V=(var(aux_SIM))/Ndelta^2;
    den_test=sqrt(V);
    KSS_test=test_num/den_test
end
%% STEP 8: AKM on FD
xy=Fdelta'*ydelta;
[psi,flag]=pcg(L,xy,1e-10,1000,pfun_);
eta=ydelta-Fdelta*psi;
dof=Ndelta-size(Fdelta,2)-1;
TSS=sum((ydelta-mean(ydelta)).^2);
R2_AKM=1-sum(eta.^2)/TSS;
adjR2_AKM=1-sum(eta.^2)/TSS*(Ndelta-1)/dof;

%% STEP 9: Graph
Tod=accumarray(match_OD,1);
Tod=Tod(match_OD);
psi_leave=mu-res./(1-Pii);
out=[ydelta,psi_leave,mu,Tod];
s=[filename 'OD_.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
end

