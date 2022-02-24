function [V, FIRST,aux_SIM2] = hutchison_trace_SE_task8(sigma_i,cluster,type_random,DO_R1,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,type_quadratic_form,K,N,J,sigma_predict,W_to_use,eps,delta)
disp('Calculating the second element of the KSS formula for standard error of a variance component in the high rank case')   
tic
NT=size(X,1)

if nargin<=16
eps=0.1;
delta=0.01;

if nargin<=17
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

    c = eps^(-2)*log(2/delta);
    VCM_TILDE	= parallel.pool.Constant((sigma_predict)); 
    sigma_i  	= sparse((1:size(sigma_i,1))',1:size(sigma_i,1),sigma_i);
    VCM_HAT 	= parallel.pool.Constant((sigma_i));    
    
    
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
               xsimul=VCM_TILDE.Value*xsimul;
               right=C_build(xsimul,X_c.Value,xx_c.Value,Lchol_c.Value,Lambda_B_c.Value,I_Lambda_P_c.Value,L_P_c.Value,type_quadratic_form,K,N,J);
               aux_SIM(s)=left'*VCM_HAT.Value*right;
    end

    left_part=2*mean(aux_SIM);
	

%Complete
    SUM=W_to_use'*sigma_predict*W_to_use; 
    V=(1/NT^2)*(4*SUM-left_part);
    FIRST=(1/NT^2)*(4*SUM); 
toc    
end
