function y = new_obj_max(x,LAMBDA)
Q_use=size(x,1)-1;

y = -1*(x(1:Q_use)'*LAMBDA*x(1:Q_use) + x(end));
