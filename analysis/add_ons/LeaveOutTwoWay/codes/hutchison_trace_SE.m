function [V, FIRST,aux_SIM2] = hutchison_trace_SE(cluster,type_random,DO_R1,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,type_quadratic_form,K,N,J,sigma_predict,W_to_use,eps,delta)
disp('Calculating the second element of the KSS formula for standard error of a variance component in the high rank case')   
tic
NT=size(X,1);

if nargin<=15
eps=0.1;
delta=0.01;

if nargin<=16
delta=0.01;
end    


%Objects for parallel world
xx_c = parallel.pool.Constant(xx);
Lchol_c =  parallel.pool.Constant(Lchol);
X_c = parallel.pool.Constant(X);
I_Lambda_P_c = parallel.pool.Constant(I_Lambda_P);
L_P_c = parallel.pool.Constant(L_P);
VCM = parallel.pool.Constant(diag(sigma_predict));
Lambda_B_c = parallel.pool.Constant(Lambda_B);


if cluster == 0
    %Simulate variance of quadratic form v'*C*v where v is multivariace normal with VCM
    %given by sigma_predict.
    NSIM=1000;
    aux_SIM=zeros(NSIM,1);
    if DO_R1==1
    aux_SIM2=zeros(NSIM,1);
    end
        parfor s=1:NSIM
               v=randn(NT,1).*sqrt(VCM.Value);
               aux=(X_c.Value)'*v;
               [coeff, flag]=pcg(xx_c.Value,aux,1e-5,1000,Lchol_c.Value,Lchol_c.Value'); 
               aux=v-X_c.Value*coeff;
               if strcmp(type_quadratic_form,'fe')    
               subtract=var(X_c.Value(:,N+1:N+J-1)*coeff(N+1:N+J-1))*(NT-1); % This is v'Bv, premultiply by NT to avoid double multiplication wrt last line.
               end
               if strcmp(type_quadratic_form,'pe')
               subtract=var(X_c.Value(:,1:N)*coeff(1:N))*(NT-1); % This is v'Bv, premultiply by NT to avoid double multiplication wrt last line.
               end
               if strcmp(type_quadratic_form,'cov')
               subtract=cov(X_c.Value(:,1:N)*coeff(1:N),X_c.Value(:,N+1:N+J-1)*coeff(N+1:N+J-1))*(NT-1); % This is v'Bv, premultiply by NT to avoid double multiplication wrt last line.
               subtract=subtract(1,2);
               end
               [aux, flag]=pcg(I_Lambda_P_c.Value,aux,1e-10,1000,L_P_c.Value,L_P_c.Value');
               aux_SIM(s)=subtract-v'*Lambda_B*aux;
               if DO_R1==1
               aux_SIM2(s) = subtract - v'*Lambda_B2*aux - lambda_1*(v'*x1bar)*(x1bar'*v);
               end 
        end
        left_part=var(aux_SIM);
        
end

%When clustering is on, we use hutchison estimator of the trace. We use the
%bounds derived in this paper: https://arxiv.org/pdf/1308.2475.pdf assuming
%that the underlying matrix is indeed PD.

if cluster == 1
    c = eps^(-2)*log(2/delta);
    VCM = parallel.pool.Constant((sigma_predict));    
    
    if strcmp(type_random,'rade')
        NSIM = floor(6*c);
    end
    
    if strcmp(type_random,'gaussian')
        NSIM = floor(8*c);
    end
   
    aux_SIM=zeros(NSIM,1);
    parfor s=1:NSIM
               if strcmp(type_random,'rade')      
                    xsimul=rand(NT,1);       
                    xsimul=1.*(xsimul>=0.5)-1.*(xsimul<0.5);
               end    
               if strcmp(type_random,'gaussian')      
                    xsimul=randn(NT,1);       
               end    
               left=C_build(xsimul,X_c.Value,xx_c.Value,Lchol_c.Value,Lambda_B_c.Value,I_Lambda_P_c.Value,L_P_c.Value,type_quadratic_form,K,N,J);
               xsimul=VCM.Value*xsimul;
               right=C_build(xsimul,X_c.Value,xx_c.Value,Lchol_c.Value,Lambda_B_c.Value,I_Lambda_P_c.Value,L_P_c.Value,type_quadratic_form,K,N,J);
               aux_SIM(s)=left'*VCM.Value*right;
    end

    left_part=2*mean(aux_SIM);
end

%Complete
    SUM=W_to_use'*sigma_predict*W_to_use; 
    V=(1/NT^2)*(4*SUM-left_part);
    FIRST=(1/NT^2)*(4*SUM); 
toc    
end
