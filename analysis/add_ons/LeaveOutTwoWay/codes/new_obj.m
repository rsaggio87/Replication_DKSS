function y = new_obj(x,LAMBDA)
Q_use=size(x,1)-1;

y = x(1:Q_use)'*LAMBDA*x(1:Q_use) + x(end);
