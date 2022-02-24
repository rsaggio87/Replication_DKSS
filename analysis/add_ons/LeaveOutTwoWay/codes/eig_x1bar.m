function [lambda,x1bar] = eig_x1bar(X,Q,lambda_eig,trace)
lambda=(lambda_eig.^2)/trace;
x1bar=X*Q;
norm=(sum(x1bar.^2)).^(0.5);
x1bar=x1bar./norm;
disp('checking, must be identity matrix')
x1bar'*x1bar
end

