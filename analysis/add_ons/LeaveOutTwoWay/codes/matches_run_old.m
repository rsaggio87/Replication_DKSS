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
out=[mu xb_leave];
s=['mat/leave_out_estimates' filename '.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
stop

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

%Now adjust denominator/numerator under the null
Bii=zeros(NT,1);
Bii_match=zeros(NT,1);
Xdesign=X;
parfor i=1:NT
mi=mm\M(i,:)';
mi=Mind*mi;
Bii(i)=var(mi)*(NT-1);
[xi flag]=pcg(xx,X(i,:)',1e-5,1000,Lchol,Lchol');
xi=Xdesign*xi;
Bii_match(i)=var(xi)*(NT-1)
end
clear Xdesign Mind
Lambda_B=spdiags(Bii,0,NT,NT);
bias_part=var(mu);
theta_match=bias_part-(1/NT)*y'*Lambda_B*eta_h;

%Now adjust numerator
[res_h, flag]=pcg(I_Lambda_P,res,1e-5,1000,L_P,L_P'); %leave-out residuals in auxiliary regression
Lambda_B=spdiags(Bii_match,0,NT,NT);
bias_part=var(X*b);
theta_XB=bias_part-(1/NT)*mu'*Lambda_B*res_h;

disp('Unbiased R2')
unbiased_R2=theta_XB/theta_match

%Inference on the numerator
A_b=X'*(X*b-mean(X*b));
aux=pcg(xx,A_b,1e-10,1000,Lchol,Lchol');
my_first_part=X*aux; %%this is B*mu
xy=Lambda_B*mu;
[ydelta_xi, flag]=pcg(I_Lambda_P,xy,1e-10,1000,L_P,L_P');
xy=X'*ydelta_xi;
[xi, flag]=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
xi=ydelta_xi-X*xi;
my_second_part=0.5*(Lambda_B*res_h+xi);
W_to_use=my_first_part-my_second_part;


%LLR
sigma_i=mu.*res_h;
Pii=diag(Lambda_P);
Bii=diag(Lambda_B);
hBest=1/NT^(1/3);
f = fit([Pii Bii],sigma_i,'lowess','Normalize','on','Span',hBest);
sigma_predict = feval(f,[Pii, Bii]);
selnan=isnan(sigma_predict);
sigma_predict(selnan)=sigma_i(selnan);

%Everything ready to make inference
NSIM=1000;
aux_SIM=zeros(NSIM,1);
   parfor s=1:NSIM
           v=randn(NT,1).*sqrt(sigma_predict);
           aux=X'*v;
           [coeff, flag]=pcg(xx,aux,1e-5,1000,Lchol,Lchol'); 
           aux=v-X*coeff;
           subtract=var(X*coeff)*(NT-1);
           [aux, flag]=pcg(I_Lambda_P,aux,1e-10,1000,L_P,L_P');
           aux_SIM(s)=subtract-v'*Lambda_B*aux;
    end
left_part=var(aux_SIM);
inner=(W_to_use).*(W_to_use).*sigma_predict;
SUM=sum(inner);
V=(1/NT^2)*(4*SUM-left_part);

%Report CI
s=['SE under q=0: ' num2str(sqrt(V))];
disp(s);
s=['CI under q=0: ' num2str((theta_XB-1.96*sqrt(V))/(theta_match)) '  '  '  ' num2str(theta_XB+1.96*sqrt(V)/(theta_match))];
disp(s);  

diary off

