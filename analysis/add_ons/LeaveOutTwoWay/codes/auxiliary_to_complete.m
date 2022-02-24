%% this is to verify calculations to compute leave out estimator. 
path(path,'~/matlab_bgl/');
%% generate Balanced Panel.
J=40;
T=5;
% very simple DGP.
%N=100;
%NT=N*T;
%id1=(1:N);
%id2=(1:N);
%id=[id1;id2];
%id=id(:);
%% very simple DGP.
N=1000;
id=(1:N)';
Tii=T;
id_py=[]
for i=1:N
%id_py=[id_py; repmat(id(i),Tii(i),1)];
id_py=[id_py; repmat(id(i),T,1)]; 
end
id=id_py;
NT=size(id,1);
firmid=1.*(id<100)+randi(J,NT,1).*(id>=100);
%DGP
NT=size(id,1);
N=max(id);
D=sparse(1:NT,id',1);
F=sparse(1:NT,firmid',1);
alpha=rand(N,1);
psi=rand(J,1);
[~,~,match_id]=unique([id firmid],'rows');
M=sparse(1:NT,match_id',1);
random_effect=randn(max(match_id),1);
y=D*alpha+F*psi+M*random_effect+randn(NT,1);
cov(D*alpha,F*psi)
year = cell2mat(accumarray(id,ones(NT,1),[],@(x){cumsum(x)}));
out=[id,firmid,year,y];
s=['src/random.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
stop
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
X=[D F*S];
xx=X'*X;
Lambda_P=sparse(NT,NT);
for i=1:max(id)
    sel=find(id==i);
    Ti=size(sel,1);
    for t=1:Ti
        aux=xx\X(sel(t),:)';    
        for tt=t:Ti
        Lambda_P(sel(tt),sel(t))=X(sel(tt),:)*aux;
        end
    end     

end
Lambda_P=Lambda_P';
Lambda_P=Lambda_P+triu(Lambda_P,1)';
I_L=speye(NT)-Lambda_P;
eta=randn(NT,1);
cav=I_L\eta;



% stop
% 
% % draw=rand(1,N);
% % firmid1=randi(J,N,1)';
% % firmid2=firmid1.*(draw<=0.95)+(randi(J,N,1))'.*(draw>0.95);
% % firmid3=firmid2.*(draw<=0.95)+(randi(J,N,1))'.*(draw>0.95);
% % firmid=[firmid1;firmid2;firmid3];
% % firmid=firmid(:);
% % size(id)
% % size(firmid)
% 
% %firmid=randi(J,NT,1);
% % %generate unbalanced Panel.
% % J=21;
% % T=6;
% % rng(1234)
% % % very simple DGP.
% % N=100;
% % id=(1:N)';
% % Tii=randi(T,N,1);
% % id_py=[]
% % for i=1:N
% % id_py=[id_py; repmat(id(i),Tii(i),1)];        
% % end
% % id=id_py;
% % NT=size(id,1);
% % 
% %%
% y=randn(NT,1);
% controls=randn(NT,1);
% prov_indicator=randn(NT,1);
% gcs = [NaN; id(1:end-1)];
% gcs = id~=gcs;
% lagfirmid=[NaN; firmid(1:end-1)];
% lagfirmid(gcs==1)=NaN; %%first obs for each worker
% [y,id,firmid,id_old,firmid_old,controls,prov_indicator] = connected_set(y,id,firmid,lagfirmid,controls,prov_indicator);
% 
% %% generate the key matrices
% NT=size(y,1);
% N=max(id);
% J=max(firmid);
% D=sparse(1:NT,id',1);
% F=sparse(1:NT,firmid',1);
% K=0;
% controls=randn(NT,K);
% JJ=F'*F;
% S=speye(J-1);
% S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
% X=[D F*S controls];
% %X=[D -F];
% xx=full(X'*X);
% S_xx_inv=xx^(-1);
% Linv=pinv(xx);
% F=F*S;
% J=size(F,2);
% alpha=randn(N,1);
% psi=rand(J,1);
% beta=[alpha;psi];
% cov(D*alpha,F*psi)
% 
% %% collapse to matches and show that Mikkel and Pat are wrong
% gcs = [NaN; id(1:end-1)];
% gcs = id~=gcs;
% lagfirmid=[NaN; firmid(1:end-1)];
% lagfirmid(gcs==1)=NaN; %%first obs for each worker
% stayer=(firmid==lagfirmid);
% stayer(gcs==1)=1;
% stayer=accumarray(id,stayer);
% T=accumarray(id,1);
% stayer=T==stayer;
% movers=stayer~=1;
% movers=movers(id);
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% A_psi=[sparse(N,N+J+K); sparse(J,N) F'*Q*F sparse(J,K); sparse(K,K+N+J)];
% P=X*S_xx_inv*X';
% B=X*S_xx_inv*A_psi*S_xx_inv*X';
% [~,~,match_id]=unique([id firmid],'rows');
% elist=index_constr(match_id,id,match_id);
% Lambda_P=diag(P)
% Lambda_B=diag(B);
% I_Lambda_P=speye(NT)-Lambda_P;
% L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
% I_Lambda_P_inv=inv(speye(NT)-Lambda_P);
% C=B-0.5*(Lambda_B*I_Lambda_P_inv*(speye(NT)-P)+(speye(NT)-P)*I_Lambda_P_inv*Lambda_B);
% Omega = rand (NT,1);
% Omega = sparse((1:NT)',(1:NT),Omega,NT,NT);
% trCCIA=trace(Omega*C*Omega*C);
% NSIM=1000;
% aux_sim=zeros(NSIM,1);
%      parfor s=1:NSIM
%                %Rademacher    
%                xsimul=rand(NT,1);
%                xsimul=1.*(xsimul>0.5)-1.*(xsimul<=0.5);
%                left=C_build(xsimul,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,'fe',K,N,J);
%                xsimul=Omega*xsimul;
%                right=C_build(xsimul,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,'fe',K,N,J);
%                aux_SIM(s)=left'*Omega*right;
%       end 
% 
% mean(aux_sim)/NT^2
% trCCIA/NT^2
% 
% [~,~,match_id]=unique([id firmid],'rows');  
% elist=index_constr(match_id,id,match_id);
% M=max(match_id);
% P=cell(M);
% for m=1:M
%     aux=X(match_id==m,:);
%     aux=xx\aux';
%     aux=X(match_id==m,:)*aux;
%     P{m}=aux;
% end
% P=P{:};
% 
% 
% 
% 
% % %%Atilde for Ftest
% % qq=5;
% % W=randn(NT,qq);
% % ww=W'*W;
% % Transform=[sparse(NT,N) F*S];
% % restrict=sparse((1:qq)',1:qq,1,qq,qq);
% % R=ww\(W'*Transform);
% % S_XX_sqrt_inv=S_xx_inv^(1/2);
% % Atilde=S_XX_sqrt_inv*(R'*R)*S_XX_sqrt_inv;
% % Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% % [Q1 lambda_eig] = eigs(Atilde,qq);
% % FunAtimesX = @(x) R'*(R*x);
% % [Q2 lambda_eig_2]=eigs(FunAtimesX,N+J-1,xx,qq);
% % lambda_eig
% % lambda_eig_2
% % W=X*Q2;
% % b=W'*y
% % sum(W(:,1).*y)
% % stop
% 
% %% Find out what is Bii for stayers
% gcs = [NaN; id(1:end-1)];
% gcs = id~=gcs;
% lagfirmid=[NaN; firmid(1:end-1)];
% lagfirmid(gcs==1)=NaN; %%first obs for each worker
% stayer=(firmid==lagfirmid);
% stayer(gcs==1)=1;
% stayer=accumarray(id,stayer);
% T=accumarray(id,1);
% stayer=T==stayer;
% movers=stayer~=1;
% movers=movers(id);
% 
% 
% %Stuff for Bii_PE
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% A_alpha=[D'*Q*D sparse(N,K+J); sparse(J+K,N+J+K)];
% Bii=zeros(NT,1);
% for i=1:NT
% xi=X(i,:);    
% Bii(i)=xi*Linv*A_alpha*Linv*xi';
% end
% quantile(Bii(~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% T=T(id);
% maxT=max(T)
% for t=1:maxT
% quantile(Bii(T==t & ~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% end 
% 
% 
% %Stuff for Bii_COV
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% A_psi_alpha=[sparse(N,N) 0.5*D'*Q*F sparse(N,K); 0.5*F'*Q*D sparse(J,J)  sparse(J,K); sparse(K,K+N+J)];
% Bii=zeros(NT,1);
% for i=1:NT
% xi=X(i,:);    
% Bii(i)=xi*Linv*A_psi_alpha*Linv*xi';
% end
% quantile(Bii(~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% for t=1:maxT
% quantile(Bii(T==t & ~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% end
% 
% 
% %Stuff for Bii_FE
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% A_psi=[sparse(N,N+J+K); sparse(J,N) F'*Q*F sparse(J,K); sparse(K,K+N+J)];
% Bii=zeros(NT,1);
% for i=1:NT
% xi=X(i,:);    
% Bii(i)=xi*Linv*A_psi*Linv*xi';
% end
% quantile(Bii(~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% for t=1:maxT
% quantile(Bii(T==t & ~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% end
% 
% %Stuff for Pii
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% Pii=zeros(NT,1);
% for i=1:NT
% xi=X(i,:);    
% Pii(i)=xi*Linv*xi';
% end
% quantile(Pii(~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% for t=1:maxT
% quantile(Pii(T==t & ~movers),[0.01 0.10 0.25 0.5 0.75 0.9 0.99])
% end
% stop
% %% Eigenvalues Calculations.
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% 
% 
% %Variance of firm effects
% A_psi=[sparse(N,N+J+K); sparse(J,N) F'*Q*F sparse(J,K); sparse(K,K+N+J)];
% S_XX_sqrt_inv=S_xx_inv^(1/2);
% Atilde=S_XX_sqrt_inv*A_psi*S_XX_sqrt_inv;
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% [Q1 lambda_eig] = eigs(Atilde,5);
% lambda_eig
% (beta'*A_psi*beta)/(NT-1)
% degree_firms=diag(F'*F);
% %degree=F'*F-(1/NT)*(degree_firms)*degree_firms';
% %A_psi_sparse=[sparse(N,N+J+K); sparse(J,N) degree sparse(J,K); sparse(K,K+N+J)];
% %[Q2 lambda_eig_2] = eigs(A_psi_sparse,xx,5);
% %%New trick of Mikkel
% Lchol=ichol(X'*X,struct('type','ict','droptol',1e-2,'diagcomp',.1));
% type_quadratic_form='fe';
% [Q3,lambda_eig_3]=eigAux(type_quadratic_form,xx,Lchol,F,D,K);
% lambda_eig_3
% 
% 
% %Re-try variance of person effects
% A_alpha=[D'*Q*D sparse(N,K+J); sparse(J+K,N+J+K)];
% Atilde=S_XX_sqrt_inv*A_alpha*S_XX_sqrt_inv;
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% [Q1 lambda_eig] = eigs(Atilde,5);
% lambda_eig
% (beta'*A_alpha*beta)/(NT-1)
% %%New trick of Mikkel
% type_quadratic_form='pe';
% [Q3,lambda_eig_3]=eigAux(type_quadratic_form,xx,Lchol,F,D,K);
% 
% lambda_eig_3
% 
% 
% %Covariance.
% A_psi_alpha=[sparse(N,N) 0.5*D'*Q*F sparse(N,K); 0.5*F'*Q*D sparse(J,J)  sparse(J,K); sparse(K,K+N+J)];
% Atilde=S_XX_sqrt_inv*A_psi_alpha*S_XX_sqrt_inv;
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% 
% [Q1 lambda_eig] = eigs(Atilde,5);
% lambda_eig
% (beta'*A_psi_alpha*beta)/(NT-1)
% type_quadratic_form='cov';
% [v,lambda_eig_3]=eigAux(type_quadratic_form,xx,Lchol,F,D,K);
% 
% lambda_eig_3
% max(max(A_psi_alpha*v-lambda_eig_3(1,1)*xx*v))
% 
% 
% 
% %For model in FD
% F=Fold;
% J=size(F,2);
% Flag=[NaN(1,J); F(1:end-1,:)];
% Flag(gcs==1,1)=NaN; %just put it for first firm, we'll remove rows where there is at least one NaN.
% Fdelta=F-Flag;
% sel=~any(isnan(Fdelta),2);
% Fdelta=Fdelta(sel,:);
% L=full(Fdelta'*Fdelta);
% 
% A=F'*Q*F;
% Linv=pinv(L);
% S_XX_sqrt_inv=Linv^(1/2);
% Atilde=S_XX_sqrt_inv*A*S_XX_sqrt_inv;
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% [Q3 lambda_eig_3] = eigs(Atilde,5);
% lambda_eig
% [Q4 lambda_eig_4] = eigs(F'*F,L,5);
% 
% 
% 
% 
% stop
% 
% 
% 
% 
% 
% 
% 
% 
% 
% S=speye(J-1);
% S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
% F=F*S;
% J=size(F,2);
% A_psi=[sparse(N,N+J+K); sparse(J,N) F'*Q*F sparse(J,K); sparse(K,K+N+J)];
% X=[D F];
% xx=X'*X;
% S_XX_sqrt=full(xx)^(1/2);
% S_XX_sqrt_inv=S_XX_sqrt^(-1);
% Atilde=S_XX_sqrt_inv*A_psi*S_XX_sqrt_inv;
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% [Q3 lambda_eig_3] = eigs(Atilde,5);
% lambda_eig_3
% 
% alpha=rand(N,1);
% psi=rand(J,1);
% psi=psi-psi(J);
% beta=randn(K,1);
% 
% X=[D,F*S,controls];
% y=D*alpha+F*S*psi(1:J-1)+controls*beta+randn(NT,1);
% gcs = [NaN; id(1:end-1)];
% gcs = id~=gcs;
% lagfirmid=[NaN; firmid(1:end-1)];
% lagfirmid(gcs==1)=NaN; %%first obs for each worker
% stayer=(firmid==lagfirmid);
% stayer(gcs==1)=1;
% stayer=accumarray(id,stayer);
% T=accumarray(id,1);
% stayer=T==stayer;
% movers=stayer~=1;
% movers=movers(id);
% %% WEIGHTING
% Dweight=F'*F;
% Q=speye(NT)-(1/NT)*ones(NT,1)*ones(1,NT);
% %% Simple checks
% xx=X'*X;
% xy=X'*y;
% tic
% Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
% toc
% b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
% bhat=b;
% ahat=b(1:N);
% ghat=b(N+1:N+J-1);
% F=F*S;
% J=size(F,2);
% %% Create the magical matrices
% A_psi=[sparse(N,N+J+K); sparse(J,N) F'*Q*F sparse(J,K); sparse(K,K+N+J)];
% A_psi_alpha=[sparse(N,N) 0.5*D'*Q*F sparse(N,K); 0.5*F'*Q*D sparse(J,J)  sparse(J,K); sparse(K,K+N+J)];
% A_alpha=[D'*Q*D sparse(N,K+J); sparse(J+K,N+J+K)];
% A_psi_2=[sparse(N,N+J+K); sparse(J,N) F'*F sparse(J,K); sparse(K,K+N+J)];
% %% Thinking about the covariance
% COV=cov(F*ghat,D*ahat);
% COV(1,2)
% (b'*A_psi_alpha*b)/(NT-1)
% A_psi_alpha_alt=[sparse(N,N) D'*Q*F sparse(N,K); sparse(J+K,J+N+K)];
% (b'*A_psi_alpha_alt*b)/(NT-1)
% %% Checking Bii now. Variance of firm effects first.
% Bii=zeros(NT,1);
% Bii_2=zeros(NT,1);
% for i=1:NT
% Xi=X(i,:);    
% Bii(i)=Xi*xx^(-1)*A_psi*xx^(-1)*Xi';
% xtilde=xx^(-1)*Xi';
% aux=xtilde(N+1:N+J);
% COV=cov(F*aux,F*aux);
% Bii_2(i)=COV(1,2)*(NT-1);
% end
% %stayers
% disp('FE')
% quantile(Bii(~movers),[0.01 0.05 0.10 0.25 0.5 0.75 0.90 0.95 0.99])
% %% Checking Bii now. Variance of person effects.
% Bii=zeros(NT,1);
% Bii_2=zeros(NT,1);
% for i=1:NT
% Xi=X(i,:);    
% Bii(i)=Xi*xx^(-1)*A_alpha*xx^(-1)*Xi';
% xtilde=xx^(-1)*Xi';
% aux=xtilde(1:N);
% COV=cov(D*aux,D*aux);
% Bii_2(i)=COV(1,2)*(NT-1);
% end
% %stayers
% disp('PE')
% quantile(Bii(~movers),[0.01 0.05 0.10 0.25 0.5 0.75 0.90 0.95 0.99])
% %% Checking Bii now. Covariance of person,worker effects.
% Bii=zeros(NT,1);
% Bii_2=zeros(NT,1);
% for i=1:NT
% Xi=X(i,:);    
% Bii(i)=Xi*xx^(-1)*A_psi_alpha*xx^(-1)*Xi';
% xtilde=xx^(-1)*Xi';
% aux=xtilde(1:N);
% COV=cov(D*aux,D*aux);
% Bii_2(i)=COV(1,2)*(NT-1);
% end
% %stayers
% disp('COV')
% quantile(Bii(~movers),[0.01 0.05 0.10 0.25 0.5 0.75 0.90 0.95 0.99])
% %% Checking Pii now.
% Bii=zeros(NT,1);
% Bii_2=zeros(NT,1);
% for i=1:NT
% Xi=X(i,:);    
% Bii(i)=Xi*xx^(-1)*Xi';
% xtilde=xx^(-1)*Xi';
% Bii_2(i)=Xi*xtilde;
% end
% disp('PII')
% quantile(Bii(~movers),[0.01 0.05 0.10 0.25 0.5 0.75 0.90 0.95 0.99])
% stop
% 
% 
% %% Now eigenvalues-eigenvectors
% degree_firms=diag(F'*F);
% degree=F'*F-(1/NT)*(degree_firms)*degree_firms';
% %degree=F'*F;
% S_XX_sqrt=full(xx)^(1/2);
% S_XX_sqrt_inv=S_XX_sqrt^(-1);
% Atilde=S_XX_sqrt_inv*A_psi*S_XX_sqrt_inv;
% Atilde=(Atilde + Atilde.')/2; %ensuring symmetry
% %lambda_1=eig(Atilde);
% %lambda_1=sort(lambda_1);
% %lambda_2=eig(full(A_psi_sparse),full(xx));
% %lambda_2=sort(lambda_2);
% %max(abs(lambda_2-lambda_1))
% A_psi_sparse=[sparse(N,N+J+K); sparse(J,N) degree sparse(J,K); sparse(K,K+N+J)];
% [Q1 lambda_eig] = eigs(Atilde,5);
% [Q2 lambda_eig_2] = eigs(A_psi_sparse,xx,5);
% lambda_eig
% lambda_eig_2
% 
% %Calculate the trace
% if 0 == 1
% trace_fe_orig=trace(Atilde*Atilde)
% Atilde=S_XX_sqrt_inv*A_alpha*S_XX_sqrt_inv;
% trace_pe_orig=trace(Atilde*Atilde)
% Atilde=S_XX_sqrt_inv*A_psi_alpha*S_XX_sqrt_inv;
% trace_cov_orig=trace(Atilde*Atilde)
% Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
% [trace_fe,trace_pe,trace_cov] = trace_Atilde_sqr(F,D,xx,Lchol)
% end
% %% Now something very silly
% v=(1/J)*[sparse(1,N) ones(1,J) sparse(1,K)]';
% %v=[sparse(1,N) 1 sparse(1,J-1) sparse(1,K)];
% A=v*v';
% Atilde=S_XX_sqrt_inv*A*S_XX_sqrt_inv;
% [Q1 lambda_eig] = eigs(Atilde,1);
% x1bar1=Q1'*S_XX_sqrt_inv*X';
% x1bar1=x1bar1';
% max(x1bar1.^2)
% stop

