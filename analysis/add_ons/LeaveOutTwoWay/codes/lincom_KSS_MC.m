function coverage_MC = lincom_KSS_MC(y,X,Z,Transform,clustering_var,Lambda_P,labels,restrict,NSIM)
%----------------------------------------------------
%%%Description
%----------------------------------------------------
% Montecarlo study for the command "lincom_KSS". See
% documentation inside "lincom_KSS" for details on inputs,
% outputs.

%----------------------------------------------------

%% READ INPUTS
got_Pii=0;
got_labels=0;
got_restrict=0;
if nargin < 9 
    NSIM=1000; %now this refers to MC draws
end
if nargin==6 && ~isempty(Lambda_P)
    got_Pii=1;
end
if nargin==7 && ~isempty(Lambda_P) && ~isempty(labels)
    got_Pii=1;
    got_labels=1;
end
if nargin==7 && isempty(Lambda_P) && ~isempty(labels)
    got_labels=1;
end
if nargin==7 && ~isempty(Lambda_P) && isempty(labels)
    got_Pii=1;
end
if nargin==8 && ~isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_labels=1;
    got_restrict=1;
end
if nargin==8 && isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_labels=1;
    got_restrict=1;
end
if nargin==8 && ~isempty(Lambda_P) && isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_restrict=1;
end
if nargin==8 && ~isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_Pii=1;
    got_labels=1;
end 
if nargin==8 && isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_labels=1;
end
if nargin==8 && ~isempty(Lambda_P) && isempty(labels) && isempty(restrict)
    got_Pii=1;
end
if nargin==9 && ~isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_labels=1;
    got_restrict=1;
end

if nargin==9 && isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict) 
    got_labels=1;
    got_restrict=1;
end
if nargin==9 && ~isempty(Lambda_P) && isempty(labels) && ~isempty(restrict) 
    got_Pii=1;
    got_restrict=1;
end
if nargin==9 && ~isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_Pii=1;
    got_labels=1;
end 
if nargin==9 && isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_labels=1;
end
if nargin==9 && ~isempty(Lambda_P) && isempty(labels) && isempty(restrict)
    got_Pii=1;
end

if nargin==10 && ~isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_labels=1;
    got_restrict=1;
end

if nargin==10 && isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict) 
    got_labels=1;
    got_restrict=1;
end
if nargin==10 && ~isempty(Lambda_P) && isempty(labels) && ~isempty(restrict) 
    got_Pii=1;
    got_restrict=1;
end
if nargin==10 && ~isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_Pii=1;
    got_labels=1;
end 
if nargin==10 && isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_labels=1;
end
if nargin==10 && ~isempty(Lambda_P) && isempty(labels) && isempty(restrict)
    got_Pii=1;
end

%% SET DIMENSIONS
n=size(X,1);
K=size(X,2);
%% Add Constant
Zold=Z;
Z=[ones(size(Z,1),1) Z];
%% PART 1: ESTIMATE HIGH DIMENSIONAL MODEL
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=X'*y;
beta=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
eta=y-X*beta;

%% PART 1B: VERIFY LEAVE OUT COMPUTATION
if got_Pii == 0 
    Lambda_P=do_Pii(X,clustering_var);
end

if got_Pii == 1 && ~isempty(clustering_var)
    nnz_1=nnz(Lambda_P);
    [~,nnz_2]=check_clustering(clustering_var);
    if nnz_1 == nnz_2
		check=['The structure of the specified Lambda_P is consistent with the level of clustering required by the user.'];
		disp(check)
	end
	if nnz_1 ~= nnz_2
		error('The user wants cluster robust inference but the Lambda_P provided by the user is not consistent with the level of clustering asked by the user. Try to omit input Lambda_P when running lincom_KSS' )
	end
end
I_Lambda_P=(speye(n,n)-Lambda_P);
eta_h=I_Lambda_P\eta; %Leave one out residual
msg = lastwarn ; 

if ~isempty(strfind(msg, 'singular'))
error('Cluster-Leave Out Residual cannot be computed. Possibly, this occurred because the level of clustering specified does not allow to compute the leave out OLS coefficient. Example: user sets cluster_var=worker_id but beta includes worker effects--> code cannot identify a worker effect after leaving out all the observations related to that worker')
end

%% PART 1C: SET UP MATRIX FOR SANDWICH FORMULA
[rows,columns]=find(Lambda_P); %read the non zero elements of Lambda_P 
aux= 0.5*(y(rows).*eta_h(columns)+y(columns).*eta_h(rows));
sigma_i=sparse(rows,columns,aux,n,n);


%% PART 2: SET UP THE TEST
r=size(Z,2);
wy=Transform*beta;
zz=Z'*Z;

%% PART 3: COMPUTE
numerator=Z\wy;
denominator=zeros(r,1);
for q=1:r
    v=sparse(q,1,1,r,1);
    v=zz\v;
    v=Z*v;
    v=Transform'*v;
    [right flag]=pcg(xx,v,1e-5,1000,Lchol,Lchol');
    left=right';    
    denominator(q)=left*(X'*sigma_i*X)*right;
end    
test_statistic=numerator./(sqrt(denominator));

%% PART 4: REPORT
if got_labels == 0 
    for q=2:r
    if q <= r    
    s=['******************************************'];
    disp(s);
    disp(s);
    disp('TRUE DGP')
    s=['******************************************'];
    disp(s); 
    s=['Linear Combination - Column Number ' num2str(q-1) ' of Z: ' num2str(numerator(q))];
    disp(s)
    s=['******************************************'];
    disp(s);
    disp(s);
    end
    end
end

if got_labels == 1 
    for q=2:r
    tell_me=labels{q-1}; 
	s=['******************************************'];
    disp(s);
    disp(s);
    disp('TRUE DGP')
    s=['******************************************'];  
    s=['******************************************'];
    disp(s); 
    s=['Linear Combination associated with '  tell_me ':  ' num2str(numerator(q))];
    disp(s)
    s=['******************************************'];
    disp(s); 
    disp(s); 
    end   
end

%% PART 5: MC

%% PART 5.1: VCM OF THE ERROR TERM
%sigma_i = abs(sigma_i); %no smoothing for the moment.
MSE = sum(eta.^2)/(n-K-1);
factor=0;
M=0;
MM=0;
if  ~isempty(clustering_var)
	[~,~,clustering_var]=unique(clustering_var);
	M=max(clustering_var);
	MM  = sparse((1:n)',clustering_var,1,n,M);
	factor = 0.1*var(y);
	random_effect = sqrt(factor)*randn(M,1);
	random_effect = MM*random_effect;
end

%% PART 5.2: SET-UP DGP
Z=Zold;
r=size(Z,2);
linear_combo_true=numerator(2:end);
linear_combo_MC=zeros(r,NSIM);
SE_linear_combo=zeros(r,NSIM);
SE_linear_combination_RES=zeros(r,NSIM);
SE_linear_combination_naive=zeros(r,NSIM);

%% PART 5.3: RUN SIMULATION
parfor s=1:NSIM
	if ~isempty(clustering_var)
		eta = (sqrt(factor)*MM*randn(M,1)) + sqrt(MSE)*randn(n,1);
	end
	if isempty(clustering_var)
		eta = sqrt(MSE)*randn(n,1);
	end
	y = X*beta+ eta;
	[~,linear_combo_MC(:,s),SE_linear_combo(:,s), SE_linear_combination_RES(:,s), SE_linear_combination_naive(:,s)] = lincom_KSS(y,X,Z,Transform,clustering_var,Lambda_P);
end
coverage_MC=zeros(r,NSIM);
oracle_coverage=zeros(r,NSIM);
coverage_MC_RES=zeros(r,NSIM);
coverage_MC_naive=zeros(r,NSIM);
for q=1:r
		aux=linear_combo_MC(q,:)';
		STD=SE_linear_combo(q,:)';
		%KSS
		UB=aux+1.96*STD;
		LB=aux-1.96*STD;
		coverage_MC(q,:)=(linear_combo_true(q)>=LB).*(linear_combo_true(q)<=UB);
		%HC1
		STD=SE_linear_combination_RES(q,:)';
		UB=aux+1.96*STD;
		LB=aux-1.96*STD;
		coverage_MC_RES(q,:)=(linear_combo_true(q)>=LB).*(linear_combo_true(q)<=UB);
		%NAIVE
		STD=SE_linear_combination_naive(q,:)';
		UB=aux+1.96*STD;
		LB=aux-1.96*STD;
		coverage_MC_naive(q,:)=(linear_combo_true(q)>=LB).*(linear_combo_true(q)<=UB);
		%ORACLE
		STD=std(aux);
		UB=aux+1.96*STD;
		LB=aux-1.96*STD;
		oracle_coverage(q,:)=((linear_combo_true(q)>=LB).*(linear_combo_true(q)<=UB));
end

%% PART 5.4: REPORT
for q=1:r
    tell_me=labels{q};
    %just focus on non-imaginary part
    sel=(imag(SE_linear_combo(q,:)) == 0); 
    SE=SE_linear_combo(q,:)';
    SE=SE(sel);
    oracle=oracle_coverage(q,:)';
    oracle=oracle(sel);
    coverage= coverage_MC(q,:);
    coverage=coverage(sel);
    s=['******************************************'];
    disp(s); 
    s=['True Value of Linear Combo associated with '  tell_me ':    ' num2str(linear_combo_true(q))];
    disp(s)
    s=['Mean Estimated Value of Linear Combo across MCs associated with '  tell_me ':    ' num2str(mean(linear_combo_MC(q,:)'))];
    disp(s)
    s=['Standard deviation of estimated Linear Combo across MCs  associated with '  tell_me ':    ' num2str(std(linear_combo_MC(q,:)'))];
    disp(s)
    s=['KSS-Avg Standard error of estimated Linear Combo across MCs associated with '  tell_me ':    ' num2str(mean(SE'))];
    disp(s)
    SE=SE_linear_combination_RES(q,:)';
    SE=SE(sel);
    s=['HC1-Avg Standard error of estimated Linear Combo across MCs associated with '  tell_me ':    ' num2str(mean(SE'))];
    disp(s)
    SE=SE_linear_combination_naive(q,:)';
    SE=SE(sel);
    s=['Naive-Avg Standard error of estimated Linear Combo across MCs associated with '  tell_me ':    ' num2str(mean(SE'))];
    disp(s)
    s=['Oracle Coverage Rate Linear Combination associated with '  tell_me ':  ' num2str(mean(oracle))];
    disp(s)
    s=['KSS-Actual Coverage Rate Linear Combination associated with '  tell_me ':  ' num2str(mean(coverage))];
    disp(s)
    coverage= coverage_MC_RES(q,:);
    coverage=coverage(sel);
	s=['HC1-Actual Coverage Rate Linear Combination associated with '  tell_me ':  ' num2str(mean(coverage))];
    disp(s)
    coverage= coverage_MC_naive(q,:);
    coverage=coverage(sel);
    s=['Naive-Actual Coverage Rate Linear Combination associated with '  tell_me ':  ' num2str(mean(coverage))];
    disp(s)
    
    s=['******************************************'];
end 

end


