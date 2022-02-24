function [y yeq,grady,gradyeq] = new_constr(x,parametri,SIGMA_inv,z_crit)

y = (parametri-x)'*SIGMA_inv*(parametri-x) -z_crit^2;

yeq=[];
grady=[];
gradyeq=[];