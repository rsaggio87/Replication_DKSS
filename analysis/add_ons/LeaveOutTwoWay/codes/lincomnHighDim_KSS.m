function test_statistic = lincomnHighDim_KSS(y,X,Z,Transform,Pii,labels)
%% Description
%This function performs inference on linear contrasts of linear
%regression coefficients parameters that is robust to many regressors and
%under heteroskedasticty. %See Remark 9 of Kline, Saggio, Solvsten (2018).

%The linear contrast is defined as v'b where v=(W'W)^(-1)W'Transform
%and b=(X'X)X'y is the OLS coefficient from the high dimensional model.
%See the example file contained in the Github repository that shows how to
%construct v in two common cases where the function lincomnHighDim_KSS can
%be particularly useful.
%% Inputs

%y: outcome from high dimensional model. Dimension n x 1.
%--
%X: regressors in the high dimensional model. Dimension n x k.
%--
%Z: matrix that defines the linear contrast. Dimension n x Q
%--
%Transform: matrix that selects and weight the object b=(X'X)X'*y.
%Dimension n x k. 
%--
%Pii: statistical leverages, i.e. X(i,:)(X'*X)^(-1)X(i,:)'. This is a NON 
%     mandatory input. If the user has already calculated this object then
%     it can be provided as input. Otherwise, the function will proceed and
%     calculate it on its own. Dimension n x 1.
%--
%labels: non-mandatory input. This is used to label the Q columns of W. It
%should be provided as a cell.
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%% Output
%test_statistic: t stat associated with the Q linear contrasts. 
%% READ INPUTS
got_Pii=0;
got_labels=0;
if nargin==5 && ~isempty(Pii)
    got_Pii=1;
end
if nargin==6 && ~isempty(Pii) && ~isempty(labels)
    got_Pii=1;
    got_labels=1;
end
if nargin==6 && isempty(Pii) && ~isempty(labels)
    got_labels=1;
end
if nargin==6 && ~isempty(Pii) && isempty(labels)
    got_Pii=1;
end
%% SET DIMENSIONS
n=size(X,1);
%% PART 1: ESTIMATE HIGH DIMENSIONAL MODEL
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=X'*y;
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
eta=y-X*b;
if got_Pii == 0 
    Pii=do_Pii(X);
end
eta_h=eta./(1-Pii); %Leave one obs out residual
sigma_i=y.*eta_h;
sigma_i=spdiags(sigma_i,0,n,n);
%% PART 2: SET UP THE TEST
Q=size(Z,2);
wy=Transform*b;
WW=Z'*Z;
%% PART 3: COMPUTE
numerator=Z\wy;
denominator=zeros(Q,1);
for q=1:Q
    v=sparse(q,1,1,Q,1);
    v=WW\v;
    v=Z*v;
    v=Transform'*v;
    [right flag]=pcg(xx,v,1e-5,1000,Lchol,Lchol');
    left=right';
    denominator(q)=left*(X'*sigma_i*X)*right;
end    
test_statistic=numerator./(sqrt(denominator));
%% PART 4: REPORT
if got_labels == 0 
    for q=1:Q
    s=['******************************************'];
    disp(s);
    disp('TESIING LINEAR RESTRICTIONS')
    s=['******************************************'];
    disp(s); 
    s=['Linear Constrast - Column Number ' num2str(q) ' of W: ' num2str(numerator(q))];
    disp(s)
    s=['Standard Error of the Linear Constract - Column Number ' num2str(q) ' of W: ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['T-stat of - for Column Number ' num2str(q) ' of W: ' num2str(test_statistic(q))];
    disp(s)
    s=['******************************************'];
    disp(s); 
    end
end
if got_labels == 1 
    for q=1:Q
    tell_me=labels{q};   
    s=['******************************************'];
    disp(s); 
    s=['Linear Constrast associated with '  tell_me ':  ' num2str(numerator(q))];
    disp(s)
    s=['SE Linear Constrast associated with '  tell_me ':  ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['T-stat: ' num2str(test_statistic(q))];
    disp(s)
    s=['******************************************'];
    disp(s); 
    end
end


