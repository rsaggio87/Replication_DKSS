function [UB,LB,C] = AM_CI(NT,lambda_1,COV_R1,bar_Beta_1,theta_2)

Q_use=size(lambda_1,1);

%Compute the Curvature of the Problem
LAMBDA								= diag(lambda_1)/NT;
parametri							= [bar_Beta_1;theta_2];
sqrt_Vb								= sqrt(COV_R1(1:Q_use,1:Q_use));
lambda_dot							= eigs(sqrt_Vb*LAMBDA*sqrt_Vb,1);
invertito							= inv(COV_R1(1:Q_use,1:Q_use));
den									= sqrt(COV_R1(Q_use+1,Q_use+1)-COV_R1(Q_use+1,1:Q_use)*invertito*COV_R1(Q_use+1,1:Q_use)');
C									= 2*abs(lambda_dot)/den;

%Search for corresponding quantile
data 								= importdata('src/AM_Tabulation5.csv'); 
data								= data.data;
Cgrid								= data(:,1);
quant_simul							= data(:,Q_use+1);
dist								= abs(C-Cgrid);
sel									= (dist==min(dist));
z_crit								= quant_simul(sel)



options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'MaxFunctionEvaluation',50000);
    
fun = @(x)new_obj(x,LAMBDA);
SIGMA_inv=inv(COV_R1);
nonlconstr = @(x)new_constr(x,parametri,SIGMA_inv,z_crit);
x0 = zeros(Q_use+1,1); % column vector

[a,LB] = fmincon(fun,x0,[],[],[],[],[],[],nonlconstr,options); 
    
fun = @(x)new_obj_max(x,LAMBDA);
[a,UB] = fmincon(fun,x0,[],[],[],[],[],[],nonlconstr,options); 

UB=-UB;    

end


