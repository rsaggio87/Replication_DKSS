function C_vec = C_build_FD(veC,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,F)
%This computes C*vec where vec is a genderal NTx1 vector and C is defined
%in KSS for a the variance of firm effects.

denom	= sum(sum(F'*F));

%Projection
xy=X'*veC;
[coeff, flag]=pcg(xx,xy,1e-5,1000,Lchol);
res_h=veC-X*coeff;
[res_h, flag]=pcg(I_Lambda_P,res_h,1e-10,1000,L_P,L_P');

eff=F*coeff;
A_b=F'*(eff-mean(eff));
  
%This completes the right part
 C_vec = construc_W_FD(veC,X,xx,Lchol,A_b,Lambda_B,I_Lambda_P,L_P,res_h);

end

