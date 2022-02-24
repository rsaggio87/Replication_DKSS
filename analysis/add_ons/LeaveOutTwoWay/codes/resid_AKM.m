function y = resid_AKM(y,id,firmid,controls)
	[~,~,id_use]=unique(id);
	[~,~,firmid_use]=unique(firmid);
	NT=size(y,1);
	D=sparse(1:NT,id_use',1); 
	F=sparse(1:NT,firmid_use',1);
	N=size(D,2);
	J=size(F,2);
	S=speye(J-1);
	S=[S;sparse(-zeros(1,J-1))];
	F=F*S;
	X=[D,F,controls]; %Design Matrix.
	xx=X'*X;
    Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
	xy=X'*y;
	b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
	y=y-X(:,N+J:end)*b(N+J:end);
end 
	