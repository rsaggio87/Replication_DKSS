function [coefficients SE]=lincom_RAKIM(y,id,firmid,lagfirmid,controls,list_lincom,sigma_i,which_firm_effect);

%Reshape sigma_i
eta_h			= sigma_i./y;
sigma_i			= (y-mean(y)).*eta_h;

%Setting up the matrices of the high-dimensional model	
NT				= size(y,1);

J				= max(firmid);
Jlag			= max(lagfirmid);
N				= max(id);

F				= sparse((1:NT)',firmid',1,NT,J);
Flag			= sparse((1:NT)',lagfirmid',1,NT,Jlag);	
D				= sparse((1:NT)',id',1,NT,N);

X				= [D F Flag controls];
K				= size(X,2);

%Setting up the matrices of the small projection model
VA              = list_lincom(:,4);
peso			= list_lincom(:,3);
peso            = sqrt(peso);
J_with_VA		= size(list_lincom,1);
sel				= list_lincom(:,2)';
if which_firm_effect == 2
tell_me 		= ['Lambda']
Transform		= [sparse(J_with_VA,N+J) 	 sparse((1:J_with_VA)',sel,peso,J_with_VA,Jlag)    sparse(J_with_VA,K-N-J-Jlag)];	
end
if which_firm_effect == 1
tell_me 		= ['Psi']
Transform		= [sparse(J_with_VA,N) 	 sparse((1:J_with_VA)',sel,peso,J_with_VA,J)    sparse(J_with_VA,K-N-J)];	
end


%Generate above and below spline
tagliami		= 3.698097;
Z				= [ones(J_with_VA,1) VA.*(VA<tagliami) VA.*(VA>=tagliami)];
Z               = peso.*Z;
r				= size(Z,2);

%%Estimate high-dimensional model
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=X'*y;
[beta]=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
eta=y-X*beta;

alpha = max(sum(abs(xx),2)./diag(xx))-2;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',alpha));
wy=(Transform*beta);
zz=Z'*Z;
coefficients=Z\wy;
denominator=zeros(r,1);
sigma_i=sparse((1:NT)',1:NT,sigma_i,NT,NT);
for q=2:r
    v=sparse(q,1,1,r,1);
    v=zz\v;
    v=Z*v;
    v=(Transform'*v);
    right=pcg(xx,v,1e-3,50000,Lchol,Lchol');
    left=right';    
    denominator(q)=left*(X'*sigma_i*X)*right;
end  
SE=sqrt(denominator);

%% PART 4: REPORT
	s=['******************************************'];
    disp(s);
    disp(s);
    disp(['INFERENCE ON LINEAR COMBINATIONS' tell_me])
    s=['******************************************'];
    disp(s);  
    disp(s); 
    for q=2:r
    if q <= r    
    s=['Linear Combination - Column Number ' num2str(q-1) ' of Z: ' num2str(coefficients(q))];
    disp(s)
    s=['Standard Error of the Linear Combination - Column Number ' num2str(q-1) ' of Z: ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['******************************************'];
    disp(s);
    end
    end
				    