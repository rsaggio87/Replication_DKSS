function [theta, V, COV_R1, gamma_sq,Fstatistic,b_1,theta_1,FIRST]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict,x1bar,lambda_1)
%This function computes the leave-out estimates of a given quadratic form
%in two way fixed effects model and associated standard error in the high
%rank case as well as the objects needed to construct a valid Confidence
%Interval in the weak identification case assuming that q=1.

%Below is a description of the inputs:

% 'type_quadratic_form': This tells the code whether we are calcualting the 
% variance of firm effects, the covariance of firm, person effects, or the 
% variance of person effects.
% 
% 'y': outcome variable.
% 
% 'J': number of firms.
% 
% 'N': number of workers.
% 
% 'X': covariates of the model (worker ids, firm ids, eventual controls).
% 
% 'xx': design matrix.
% 
% 'Lchol': chol of xx.
% 
% 'bias_part': plugin estimate of the variance decomposition parameter.
% 
% 'eta_h': leave-out residuals.
% 
% 'I_Lambda_P': I-Lambda_P where I is an identity matrix.
% 
% 'L_P': chol of I_Lambda_P.
% 
% 'Lambda_B': diag(Bii).
% 
% 'W_to_use': This is C*y=B*y-0.5(Lambda_B*eta_h+xi_hat) where xi_hat is 
%the residual when regressing X onto (1-Lambda_p)^(-1)*Lambda_B*y. The
%matrix C is defined below equation (9) in KSS. This object was computed in
%the auxiliary function 'construc_W'.

%'my_first_part:' This is B*y=B*Y=X(X'X)^(-1)A*(X'X)^(-1)X'Y This object was computed in
%the auxiliary function 'construc_W'.

%'sigma_predict': \tilde{\sigma}_{i}^2, see section on high rank inference.

%'x1bar': what is denoted as \mathsf{w}_{i1} in KSS.

%'lambda_1': largest eigenvalue of the matrix called Atilde in KSS. This is
%not divided by the total number of observations. 

%Below is a description of the outputs:

%'theta': leave out estimate of the associated variance decomposition
%parameter

%'V': sampling variance of theta.

%'COV_R1': \{Sigma}_1 as defined in section 5.3

%'b_1': b_{1} as defined in section 5.3

%'theta_1': \theta_{1} as defined in section 5.3

%See end of code for definition of 'gamma_sq' and 'Fstatistic'.
K = size(X,2)-N-J+1;
do_SE=1;
DO_R1=0;
V=0;
COV_R1=0;
gamma_sq=0;
Fstatistic=0;
b_1=0;
theta_1=0;
FIRST=0;


if nargout==1
    do_SE=0;
end

if nargout>2 && nargin >16
DO_R1=1;  
aux=lambda_1*(x1bar.^2);
aux=spdiags(aux,0,size(eta_h,1),size(eta_h,1));
Lambda_B2=Lambda_B-aux;
end
if nargout<=2 || nargin <=16
x1bar=0;
lambda_1=0;
DO_R1=0;
Lambda_B2=sparse(size(eta_h,1),size(eta_h,1));    
end   

%% Set up dimensions
NT=size(eta_h,1);
%% STEP 1: Leave out estimate of the quadratic form
theta=bias_part-(1/NT)*y'*Lambda_B*eta_h;

%% Step 2: Construct SE in the high rank case
%Here we focus on estimation of the sampling variance of theta in the high
%rank case, see end of page 25. We first focus on computing the right part
%of the formula. That part corresponds to the variance of a quadratic form
%with a vector of normal random variables each with VCM variance 
%\tilde{\Omega}=\sigma_predict. The function "hutchison_trace_SE" extends
%this idea to the case when \tilde{\Omega} accounts for the possibility 
%of serial correlation within match.

%In the heteroskedatic case, we approximate this object via simulations 
%We also proceed in a similar way to compute the object needed to compute
%inference in the case in which q=1.
if do_SE == 1
    cluster=1-isdiag(sigma_predict);
    type_random='rade'; %%this will only be read if cluster==1
    
    if DO_R1 == 0
    [V, FIRST] = hutchison_trace_SE(cluster,type_random,DO_R1,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,type_quadratic_form,K,N,J,sigma_predict,W_to_use);
    end
    
    if DO_R1 == 1
    [V,FIRST,aux_SIM2] = hutchison_trace_SE(cluster,type_random,DO_R1,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,type_quadratic_form,K,N,J,sigma_predict,W_to_use);
    end
    
end
    %% Step 3: Standard Errors for case when q=1
if DO_R1==1 && cluster == 0
    %First thing I need is to form C_2*y=B_2*y-0.5(Lambda_B2*eta_h+xi_hat_2)
    %where:  
    %xi_hat_2 is the residual when regressing X onto
    %(1-Lambda_p)^(-1)*Lambda_B2*y

    %B_2=B-lambda_1*x1bar*x1bar'

    %Lambda_B2=Lambda_B-diag(lambda_1*x1bar*x1bar') --> defined at the
    %beginning.



    %start with the first piece: 
    %B_2*y=(B-lambda_1*x1bar*x1bar')y

    %Recall that B*y=my_first_part
    b_1=sum(x1bar.*y);
    first_piece=my_first_part-(lambda_1*x1bar)*b_1;


    %Now the second bit
    xy=Lambda_B2*y;
    [ydelta_xi, flag]=pcg(I_Lambda_P,xy,1e-10,1000,L_P,L_P');
    xy=X'*ydelta_xi;
    [b, flag]=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
    %                 if flag > 0
    %                     error(['PCG FLAG: ' num2str(flag)])
    %                 end
    xi=ydelta_xi-X*b;
    second_piece=0.5*(Lambda_B2*eta_h+xi);
    W_2=first_piece-second_piece; %this gives me W_2(i)=\sum_{l \neq i}C_{il}y_{l}

    %Let's set-up things to calculate the VCM
    COV_R1=zeros(2,2);

    %Start with COV(1,1)
    %COV_R1(1,1)=sum(x1bar.^2.*sigma_i); %do not divide
    COV_R1(1,1)=sum(x1bar.^2.*sigma_predict); %do not divide

    %Next is COV(1,2)
    COV_R1(1,2)=2*sum(x1bar.*sigma_predict.*W_2)/(NT);
    COV_R1(2,1)=COV_R1(1,2);

    %More complicated step is COV(2,2). Similar structure of calculations as
    %performed when q=0 --> variance of a quadratic form with mean zero normal
    %vector. See the calculations computed there (inefficient to write another
    %parfor).

    second_bit=var(aux_SIM2);
    COV_R1(2,2)=(4*sum(W_2.*W_2.*sigma_predict)-second_bit)/(NT^2);

    %Gamma Squared
    gamma_sq=((lambda_1^2/NT^2)*(COV_R1(1,1)^2))/(COV_R1(2,2));

    %Fstatistic
    Fstatistic=(b_1^2)/COV_R1(1,1);

    %For AM Method
    theta_1=theta-(lambda_1/(NT))*(b_1^2-COV_R1(1,1));
end
end    


