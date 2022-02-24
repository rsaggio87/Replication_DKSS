function [Pii, Bii] = eff_res_FD_shrink(elist,Fdelta,L,pfun_)
%Computation of Pii, Bii for variance of firm effects in two-way fixed
%effects model in a model that partials out worker effects by taking first
%differences.
    
%% Auxiliary Step
    J = size(Fdelta,2);%number of nodes in the graph
    M= size(elist,1);
    
%Preallocate    
    output = zeros(M,1);
    output_B=zeros(M,1);  
    
%% Finding leverages and cross products using naive but exact parfor.
    tic  
        parfor m=1:M
            Bleft=sparse(elist(m,1:2),1,[-1 1],J,1);
            Bright=sparse(elist(m,3:4),1,[-1 1],J,1);
            [xtilde flag]=pcg(L,Bright,1e-05,300,pfun_);
            %Lambda_P
            output(m)=Bleft'*xtilde;
            %Lambda_B
            [xtilde_left flag]=pcg(L,Bleft,1e-05,300,pfun_);
            aux=cov(xtilde_left,xtilde)*(J-1);
            output_B(m)=aux(1,2);
        end
    disp('Time to Compute Leverage Scores (Exact)')
    toc
%% Export
Pii=output;
Bii=output_B;
end 
   
    
