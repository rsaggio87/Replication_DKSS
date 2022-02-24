%Routine to verify calculations. 

load([filename '_completed']) 


%verify the leave out formulas
J=J-1;
F=F*S;
n_matches=max(clustering_var);
try_out=100;
max_diff=zeros(try_out,1);


parfor s=1:try_out
    i=randi(n_matches,1,1);
    sel=find(clustering_var==i);
    Ti=size(sel,1);
    for t=1:Ti
        for tt=t:Ti
        aux=xx\X(sel(t),:)';
        Pii=X(sel(tt),:)*aux;
        pe=X(:,1:N)*aux(1:N,:);
        fe=X(:,N+1:end)*aux(N+1:end,:);
        COV=cov(pe,fe)*(NT-1);
        aux1=max(max(abs(((Pii)-Lambda_P(sel(t),sel(t))))));
        aux2=max(max(abs((COV(2,2))-Lambda_B_fe(sel(t),sel(t)))));
        if n_of_parameters == 3
        aux4=max(max(abs((COV(1,1)-Lambda_B_pe(sel(t),sel(t))))));
        end
        if n_of_parameters >= 2
        aux3=max(max(abs((COV(1,2)-Lambda_B_cov(sel(t),sel(t))))));
        end
        if n_of_parameters == 3
        max_diff_aux=max([aux1;aux2;aux3;aux4]);
        end
        if n_of_parameters == 2
        max_diff_aux=max([aux1;aux2;aux3]);
        end
        if n_of_parameters == 1
        max_diff_aux=max([aux1;aux2]);
        end
        max_diff(s)=max_diff_aux
        end
    end
end    
max(max_diff)


%verify computation of the standard error (only for small datasets).
if 1 == 1
if do_SE == 1
xx=full(X'*X);
S_xx_inv=xx^(-1);
P=X*S_xx_inv*X';
I_Lambda_P=speye(NT)-Lambda_P;
I_Lambda_P_inv=inv(speye(NT)-Lambda_P);
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);


Lambda_B=Lambda_B_fe;
A_psi=[sparse(N,N+J+K); sparse(J,N) F'*Q*F sparse(J,K); sparse(K,K+N+J)];
B=X*S_xx_inv*A_psi*S_xx_inv*X';
C=B-0.5*(Lambda_B*I_Lambda_P_inv*(speye(NT)-P)+(speye(NT)-P)*I_Lambda_P_inv*Lambda_B);
PSI=(y'*C*y)/NT
traccia=2*trace(C*sigma_predict*C*sigma_predict);
TRACCIA=traccia/NT^2
W_to_use=C*y; 
SUM=W_to_use'*sigma_predict*W_to_use;
SOMMA=SUM/NT^2
V=(4*SUM-traccia)/NT^2;
sqrt(V)
SE_sigma2_psi
stop

if n_of_parameters >= 2
Lambda_B=Lambda_B_cov;
A_psi_alpha=[sparse(N,N) 0.5*D'*Q*F sparse(N,K); 0.5*F'*Q*D sparse(J,J)  sparse(J,K); sparse(K,K+N+J)];
B=X*S_xx_inv*A_psi_alpha*S_xx_inv*X';
C=B-0.5*(Lambda_B*I_Lambda_P_inv*(speye(NT)-P)+(speye(NT)-P)*I_Lambda_P_inv*Lambda_B);
traccia=2*trace(C*sigma_predict*C*sigma_predict);
W_to_use=C*y; 
SUM=W_to_use'*sigma_predict*W_to_use;
V=(4*SUM-traccia)/NT^2;
sqrt(V)
SE_sigma_psi_alpha
end

if n_of_parameters == 3
A_alpha=[D'*Q*D sparse(N,K+J); sparse(J+K,N+J+K)];
Lambda_B=Lambda_B_pe;
B=X*S_xx_inv*A_alpha*S_xx_inv*X';
C=B-0.5*(Lambda_B*I_Lambda_P_inv*(speye(NT)-P)+(speye(NT)-P)*I_Lambda_P_inv*Lambda_B);
traccia=2*trace(C*sigma_predict*C*sigma_predict);
W_to_use=C*y; 
SUM=W_to_use'*sigma_predict*W_to_use;
V=(4*SUM-traccia)/NT^2;
sqrt(V)
SE_sigma2_alpha
end
end
end
