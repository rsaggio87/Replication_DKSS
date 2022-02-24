clc
clear
s=['mat/COMPLETE_Belluno'];
load(s)
%%%%%%%%%%%%%LOG FILE / Options
logname=['../logs/leave_one_out_complete' filename '_Idea1match_model.log'];
system(['rm ' logname])
diary(logname)

%Create leave out estimates
Pii=diag(Lambda_P);
xb_leave=xb-Pii.*eta_h;


%Start by storing the sigma_i
I_Lambda_P=(speye(NT,NT)-Lambda_P);
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));

%create match indicators.
[~,~,match_id]=unique([id firmid],'rows');
M=sparse(1:NT,match_id',1);
m=max(match_id); 
mm=M'*M;
Mind=M;

%Obtain matched means.
mu=M\y;
mu=M*mu;

%Now report naive statistics
xy=X'*mu;
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
res=mu-X*b;
dof=NT-size(X,2)-1;
TSS=sum((mu-mean(mu)).^2);
R2=1-sum(res.^2)/TSS;
adjR2=1-sum(res.^2)/TSS*(NT-1)/dof;
disp('R2 when regressing matched effects on worker, firm dummies')
R2
disp('AdjR2 when regressing matched effects on worker, firm dummies')
adjR2

%Now KSS the R2
Bii=zeros(NT,1);
Xdesign=X;
parfor i=1:NT
[xi flag]=pcg(xx,X(i,:)',1e-5,1000,Lchol,Lchol');
xi=Xdesign*xi;
Bii(i)=var(xi)*(NT-1)
end
clear Xdesign
Lambda_B=spdiags(Bii,0,NT,NT);
bias_part=var(xb);
theta_match=bias_part-(1/NT)*y'*Lambda_B*eta_h;

%Out to Stata
out=[mu xb_leave repmat(theta_match,NT,1)];
s=['mat/leave_out_estimates' filename '.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 

diary off

