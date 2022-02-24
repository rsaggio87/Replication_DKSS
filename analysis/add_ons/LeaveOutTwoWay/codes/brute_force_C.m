function [C Cq] = brute_force_C(X,F,lambda_1,x1bar)



n=size(X,1);

S_xx = X'*X;
AUX =  S_xx\X';
P 	 = X*AUX;
B	 = F*AUX;
B	 = (B-mean(B));
B	 = F'*B;
B	 = AUX'*B;
clear AUX

Lambda_B = diag(diag(B));
Lambda_P = diag(diag(P));
M		 = speye(n)-P;
I_P_inv	 = diag(M);
I_P_inv	 = 1./I_P_inv;
I_P_inv	 = sparse((1:n),1:n,I_P_inv,n,n);

C 		 		= B - 0.5*(Lambda_B*I_P_inv*M+M*I_P_inv*Lambda_B);


Lambda_B2		= lambda_1*(x1bar.^2);
Lambda_B2		= spdiags(Lambda_B2,0,size(X,1),size(X,1));
Lambda_B2		= Lambda_B-Lambda_B2;

Cq       = B-lambda_1*x1bar*x1bar'- 0.5*(Lambda_B2*I_P_inv*M+M*I_P_inv*Lambda_B2);

end

    