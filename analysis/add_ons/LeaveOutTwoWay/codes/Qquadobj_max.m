function [y,grady] = Qquadobj_max(x,Q,f,c)
y = -1.*(1/2*x'*Q*x + f'*x + c);
if nargout > 1
    grady = -Q*x + f;
end
end