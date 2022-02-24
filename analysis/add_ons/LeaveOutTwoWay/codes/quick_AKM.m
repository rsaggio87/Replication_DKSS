function fe = quick_AKM(y,id,firmid,id_orig,firmid_orig,filename)
	NT=size(y,1);
	D=sparse(1:NT,id',1); 
	F=sparse(1:NT,firmid',1);
	N=size(D,2);
	J=size(F,2);
	S=speye(J-1);
	S=[S;sparse(-zeros(1,J-1))];
	F=F*S;
	X=[D,F]; %Design Matrix.
	xx=X'*X;
    Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
	xy=X'*y;
	b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
	ahat=b(1:N);
	ghat=b(N+1:N+J-1);
	fe=F*ghat;
end 
	