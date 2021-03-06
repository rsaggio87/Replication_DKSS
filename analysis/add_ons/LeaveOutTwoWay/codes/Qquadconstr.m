function [y,yeq,grady,gradyeq] = Qquadconstr(x,H,k,d,parametri)
jj = length(H); % jj is the number of inequality constraints
y = zeros(1,jj);
for i = 1:jj
    y(i) = 1/2*(parametri-x)'*H{i}*(parametri-x) + k{i}'*x + d{i};
end
yeq = [];
    
if nargout > 2
    grady = zeros(length(x),jj);
    for i = 1:jj
        grady(:,i) = -H{i}*x + k{i};
    end
end
gradyeq = [];