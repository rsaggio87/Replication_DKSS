load([filename '_completed']) 

Pii=zeros(NT,1);
Bii=zeros(NT,1);
parfor i=1:NT
aux=xx\X(i,:)'
Pii(i)=X(i,:)*aux;
fe=X(:,N+1:end)*aux(N+1:end);
Bii(i)=var(fe)*(NT-1);
end    
xy=X'*y;
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');   
xb=X*b;
eta=y-xb;
eta_h=eta./(1-Pii); %Leave one out residual
sigma_i=y.*eta_h;
fe=F*S*ghat;

sigma2_psi
estimator=var(fe)-(1/NT)*sum(Bii.*sigma_i)