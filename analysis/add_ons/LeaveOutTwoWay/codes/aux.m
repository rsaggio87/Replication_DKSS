clc
clear
for ss=1:3
s=['mat/after_step_3_completeBelluno'];
load(s)
%%%%%%%%%%%%%LOG FILE / Options
subsample_llr_fit=ss;
logname=['../logs/leave_one_out_complete' filename '_Method' num2str(subsample_llr_fit) '.log'];
system(['rm ' logname])
diary(logname)

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
filename
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%% STEP 4: DIAGNOSTICS
%This part is only computed provided that the option 'eigen_diagno' is
%turned on. In this part, we calculate the squared eigenvalue ratio and
%Lindeberg conditions that constitute the key conditions to verify the
%validity of Theorem 1. Lindeberg condition assumes q=1.

%Future releases of this function will allow for faster computation of the 
%eigenvalue conditions. In particular, the step where we build and store 
%the A matrix for calculation of the eigenvectors can be sidestepped.

if eigen_diagno==1 
EIG_NORM=zeros(3,n_of_parameters);
max_x1bar_sq=zeros(n_of_parameters,1);
lambda_1=zeros(n_of_parameters,1);
SUM_EIG=zeros(n_of_parameters,1);
x1bar_all=zeros(NT,n_of_parameters);
tic    
%Begin by calcuting the sum of the squared of the eigenvalues for the corresponding Atilde
disp('Calculating Sum of Squared Eigenvalues via Hutchinson...')

if n_of_parameters==1
    [trace_fe]=trace_Atilde_sqr(X,F*S,D,xx,Lchol);
end

if n_of_parameters==2
    [trace_fe, trace_cov]=trace_Atilde_sqr(X,F*S,D,xx,Lchol);
end

if n_of_parameters==3
    [trace_fe,trace_cov,trace_pe]=trace_Atilde_sqr(X,F*S,D,xx,Lchol);
end

if n_of_parameters==1
    SUM_EIG(1)=trace_fe;
end

if n_of_parameters==2
    SUM_EIG(1)=trace_fe;
    SUM_EIG(2)=trace_pe;
end

if n_of_parameters==3
    SUM_EIG(1)=trace_fe;
    SUM_EIG(2)=trace_pe;
    SUM_EIG(3)=trace_cov;
end
disp('done') 

%Issue: we use a trick in order to compute the eigenvalues/vectors
%associated with the matrix Atilde without having to store the matrix
%Atilde in memory

for pp=1:n_of_parameters
  
    if pp == 1
    type_quadratic_form='fe';
    entry=['Calculating Diagnostic for Variance of Firm Effects...'];
    end
    
    if pp == 2
    type_quadratic_form='cov';
    entry=['Calculating Diagnostic for Covariance of Firm, Person Effects...'];
    end
    
    if pp == 3
    type_quadratic_form='pe';
    entry=['Calculating Diagnostic for Variance of Person Effects...'];
    end
    
    %Calculate
    disp(entry)
    [Q,lambda_eig] = eigAux(type_quadratic_form,xx,Lchol,F*S,D,K);
    lambda_eig=diag(lambda_eig);
    lambda_1(pp)=lambda_eig(1); 
    [EIG_NORM(:,pp),x1bar_all(:,pp)] = eig_x1bar(X,Q,lambda_eig,SUM_EIG(pp)); %Note: What we label as x1bar in the code corresponds to w_i1 in KSS.
    max_x1bar_sq(1)=max(x1bar_all(:,pp));
    
    
end
disp('checking, must report 1')
sum(x1bar_all.^2,1)  
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
disp('Running AKM...')
tic
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
R2=1-sum(eta.^2)/TSS;
adjR2=1-sum(eta.^2)/TSS*(NT-1)/dof;
COV=cov(fe,pe);
sigma_2_psi_AKM=COV(1,1);
sigma_2_alpha_AKM=COV(2,2);
sigma_alpha_psi_AKM=COV(1,2);

%% STEP 5.2: Leave One Out Estimation.
I_Lambda_P=(speye(NT,NT)-Lambda_P);
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[eta_h, flag]=pcg(I_Lambda_P,eta,1e-5,1000,L_P,L_P'); %leave-out residuals

%Preallocate vectors with results
theta=zeros(n_of_parameters,1);
V_theta=zeros(n_of_parameters,1);
COV_R1=zeros(2,2,n_of_parameters);
gamma_sq=zeros(n_of_parameters,1);
F_stat=zeros(n_of_parameters,1);
b_1=zeros(n_of_parameters,1);
theta_1=zeros(n_of_parameters,1);

%Loop over parameters to be estimated
for pp=1:n_of_parameters
    
    if eigen_diagno == 1
        x1_bar=x1bar_all(:,pp);
    end
    
    if pp == 1 %Variance of Firm Effects
        type_quadratic_form='fe';
        Lambda_B=Lambda_B_fe;
        bias_part=sigma_2_psi_AKM; 
        if no_controls == 0
            A_b=[zeros(N,1); S'*F'*(fe-mean(fe)); zeros(K,1)];    
        end 
        if no_controls == 1
            A_b=[zeros(N,1); S'*F'*(fe-mean(fe))];
        end
        entry=['Calculating Leave out, Variance of Firm Effects...'];
    end


    if pp == 2 %CoVariance of Person,Firm Effects
        type_quadratic_form='cov';
        Lambda_B=Lambda_B_cov;
        bias_part=sigma_alpha_psi_AKM; 
        if no_controls==0
        A_b=[0.5*D'*(fe-mean(fe)); 0.5*S'*F'*(pe-mean(pe)); zeros(K,1)];
        end
        if no_controls==1
        A_b=[0.5*D'*(fe-mean(fe)); 0.5*S'*F'*(pe-mean(pe))];
        end
        entry=['Calculating Leave out, Covariance of Person, Firm Effects...'];
    end

    if pp == 3 %Variance of Person Effects
        type_quadratic_form='pe';
        Lambda_B=Lambda_B_pe;
        bias_part=sigma_2_alpha_AKM; 
        if no_controls==0
            A_b=[D'*(pe-mean(pe)); zeros(J-1,1); zeros(K,1)];
        end
        if no_controls==1
            A_b=[D'*(pe-mean(pe)); zeros(J-1,1)];
        end
        entry=['Calculating Leave out, Variance of Person Effects...'];
    end

%Tell me what are we doing
     disp(entry)
   
%Auxiliary
    if do_SE  == 1 
        [W_to_use, my_first_part] = construc_W(y,X,xx,Lchol,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);
    end
%Compute the non-parametric fit needed for standard error estimation in 
%high rank case.
    if do_SE == 1  
        sigma_predict= llr_fit(Lambda_P,Lambda_B,y,eta_h,subsample_llr_fit,K,movers,T);
    end

%Run
    if eigen_diagno == 1  && do_SE == 0
        [theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B);
    end

    if eigen_diagno == 0  && do_SE == 0
        [theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B);
    end
    
    if eigen_diagno == 0  && do_SE == 1
        [theta(pp), V_theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict);
    end
    if eigen_diagno == 1  && do_SE == 1
        [theta(pp), V_theta(pp), COV_R1(:,:,pp), gamma_sq(pp), F_stat(pp), b_1(pp), theta_1(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict,x1_bar,lambda_1(pp));
    end
end

clear A_b Lambda_B W_to_use I_Lambda_P L_P

if n_of_parameters==1 && do_SE==0
    sigma2_psi=theta(1);
end

if n_of_parameters==2 && do_SE==0
    sigma2_psi=theta(1);
    sigma_psi_alpha=theta(2);
end   

if n_of_parameters==3 && do_SE==0
    sigma2_psi=theta(1);
    sigma_psi_alpha=theta(2);
    sigma2_alpha=theta(3);
end


if n_of_parameters==1 && do_SE==1
    sigma2_psi=theta(1);
    SE_sigma2_psi=sqrt(V_theta(1));
end

if n_of_parameters==2 && do_SE==1
    sigma2_psi=theta(1);
    SE_sigma2_psi=sqrt(V_theta(1));
    sigma_psi_alpha=theta(2);
    SE_sigma_psi_alpha=sqrt(V_theta(2));
end   

if n_of_parameters==3 && do_SE==1
    sigma2_psi=theta(1);
    SE_sigma2_psi=sqrt(V_theta(1));
    sigma_psi_alpha=theta(2);
    SE_sigma_psi_alpha=sqrt(V_theta(2));
    sigma2_alpha=theta(3);
    SE_sigma2_alpha=sqrt(V_theta(3));
end  

%R2
if n_of_parameters==3
    explained_var_leave_out=(sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var(y);
end




%% %% STEP 5.3: Andrews
if andrews_estimates==1 && no_controls==0    
[var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D,controls);
end
if andrews_estimates==1 && no_controls==1    
[var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D);
end


%% STEP 6: REPORTING
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
s=['Maximum Leverage: ' num2str(max(diag(Lambda_P)))];
disp(s);
s=['-*-*-*-*-*-*AKM*-*-*-*'];
disp(s)
s=['Variance of Firm Effects: ' num2str(sigma_2_psi_AKM)];
disp(s)
s=['Covariance of Firm and Person Effects: ' num2str(sigma_alpha_psi_AKM)];
disp(s);
s=['Variance of Person Effects: ' num2str(sigma_2_alpha_AKM)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str((corr(fe,pe)))];
disp(s);
s=['R2: ' num2str(R2)];
disp(s);
if andrews_estimates == 1
s=['-*-*-*-*-*-*ANDREWS*-*-*-*'];
disp(s);
s=['Variance of Firm Effects: ' num2str(var_corrected_fe)];
disp(s);
s=['Covariance of Firm and Person Effects: ' num2str(var_corrected_cov)];
disp(s);
s=['Variance of Person Effects: ' num2str(var_corrected_pe)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str(var_corrected_cov/(sqrt(var_corrected_fe)*sqrt(var_corrected_pe)))];
disp(s);
s=['Total Explained Variation: ' num2str(adjR2)];
disp(s);
end
s=['-*-*-*-*-*-*LEAVE ONE OUT*-*-*-*'];
disp(s)
s=['Variance of Firm Effects: ' num2str(sigma2_psi)];
disp(s);
if n_of_parameters>=2
s=['Covariance of Firm and Person Effects: ' num2str(sigma_psi_alpha)];
disp(s);
end
if n_of_parameters==3
s=['Variance of Person Effects: ' num2str(sigma2_alpha)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str(sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)))];
disp(s);
s=['Total Explained Variation: ' num2str(explained_var_leave_out)];
disp(s);
end
%% Alternative Way to Show Output
s=['-*-*-*-*-*-*Variance of Firm Effects-*-*-*-*-*-*'];
disp(s)
s=['AKM: ' num2str(sigma_2_psi_AKM)];
disp(s)
if andrews_estimates == 1
s=['Andrews: ' num2str(var_corrected_fe)];
disp(s)
end
s=['Leave One Out: ' num2str(sigma2_psi)];
disp(s)
if do_SE == 1
s=['Leave One Out SE: ' num2str(SE_sigma2_psi)];
disp(s)
end
if n_of_parameters>=2
    s=['-*-*-*-*-*-*Covariance of Firm,Person Effects-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: ' num2str(sigma_alpha_psi_AKM)];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(var_corrected_cov)];
    disp(s);
    end
    s=['Leave One Out: ' num2str(sigma_psi_alpha)];
    disp(s);
    if do_SE==1
    s=['Leave One Out SE: ' num2str(SE_sigma_psi_alpha)];
    disp(s)
    end
end
if n_of_parameters>=3
    s=['-*-*-*-*-*-*Variance of Person Effects-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: ' num2str(sigma_2_alpha_AKM)];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(var_corrected_pe)];
    disp(s);
    end
    s=['Leave One Out: ' num2str(sigma2_alpha)];
    disp(s);
    if do_SE==1
    s=['Leave One Out SE: ' num2str(SE_sigma2_alpha)];
    disp(s)
    end
    s=['-*-*-*-*-*-*Correlation of Person,Firm Effects-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: ' num2str((corr(fe,pe)))];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(var_corrected_cov/(sqrt(var_corrected_fe)*sqrt(var_corrected_pe)))];
    disp(s);
    end
    s=['Leave One Out: ' num2str(sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)))];
    disp(s);
    s=['-*-*-*-*-*-*R2-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: '  num2str(R2)];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(adjR2)];
    disp(s);
    end
    s=['Leave One Out: ' num2str(explained_var_leave_out)];
    disp(s);
end
%% Focus on Diagnostics
%Note: The eigenvalue ratios (and sum of squared eigenvalues) 
%      reported in Table 3 of KSS have been calculated exactly whereas here
%      we report the results obtained via simulations. 
%      The differences are neglible as one can see from below.
if eigen_diagno==1
    for pp=1:n_of_parameters
        
        if pp == 1
            title='Diagnostics on Variance of Firm Effects';
        end  
        if pp == 2
            title='Diagnostics on CoVariance of Person, Firm Effects';
        end
        if pp == 3
            title='Diagnostics on Variance of Person Effects';
        end

        s=['*********************' title '*********************'];
        disp(s);
        s=['ratio of eigenvalues: '];
        disp(s)
        EIG_NORM(1:3,pp)
        s=['Weak Lindeberg Condition: ' num2str(max_x1bar_sq(pp))];
        disp(s)
        s=['Sum of squared eigenvalues: ' num2str(SUM_EIG(pp)/NT^2)];
        disp(s)
        s=['Variance of b_1:  ' num2str(COV_R1(1,1,pp))];
        disp(s);
        s=['Variance of \hat{\theta}_1: ' num2str(COV_R1(2,2,pp))];
        disp(s);
        s=['Covariance of (b_1,\theta_1): ' num2str(COV_R1(1,2,pp))];
        disp(s);
        s=['Correlation of (b_1,\theta_1): ' num2str((COV_R1(1,2,pp))/(sqrt(COV_R1(2,2,pp))*sqrt(COV_R1(1,1,pp))))];
        disp(s);
        s=['gamma squared: ' num2str(gamma_sq(pp))];
        disp(s);
        s=['Fstatistic: ' num2str(F_stat(pp))];
        disp(s);
        s=['******************************************'];
        disp(s);
    end
end
%% Focus on Inference
if do_SE==1
    for pp=1:n_of_parameters
        
        if pp == 1
            title='Inference on Variance of Firm Effects';
        end
        
        if pp == 2
            title='Inference on CoVariance of Person, Firm Effects';
        end
        
        if pp == 3
            title='Inference on Variance of Person Effects';
        end
        
        s=['*********************' title '*********************'];
        disp(s);
        s=['SE under q=0: ' num2str(sqrt(V_theta(pp)))];
        disp(s);
        s=['CI under q=0: ' num2str(theta(pp)-1.96*sqrt(V_theta(pp))) '  '  '  ' num2str(theta(pp)+1.96*sqrt(V_theta(pp)))];
        disp(s);   
        
        if eigen_diagno==1
                [UB,LB,C]=AM_CI(NT,lambda_1(pp),gamma_sq(pp),COV_R1(:,:,pp),b_1(pp),theta_1(pp));
                s=['CI under q=1: ' num2str(LB) '  '  '  ' num2str(UB)];
                disp(s);  
                s=['Curvature: ' num2str(C)];
                disp(s); 
                s=['******************************************'];
                disp(s); 
        end
    end
end
diary off
end