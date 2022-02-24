function C_vec = C_build(veC,X,xx,Lchol,Lambda_B,I_Lambda_P,L_P,type_quadratic_form,K,N,J)
%This computes C*vec where vec is a genderal NTx1 vector and C is defined
%in KSS for a given type of quadratic form specified in type_quadratic_form

%Projection
xy=X'*veC;
[coeff, flag]=pcg(xx,xy,1e-5,1000,Lchol,Lchol');
res_h=veC-X*coeff;
[res_h, flag]=pcg(I_Lambda_P,res_h,1e-10,1000,L_P,L_P');


%Auxiliary step
if strcmp(type_quadratic_form,'fe') 
   eff=X(:,N+1:N+J-1)*coeff(N+1:N+J-1,1);
   
   if K>0
      A_b=[zeros(N,1); X(:,N+1:N+J-1)'*(eff-mean(eff)); zeros(K,1)];
   end
   
   if K==0
      A_b=[zeros(N,1); X(:,N+1:N+J-1)'*(eff-mean(eff))];
   end    
               
end
               
  if strcmp(type_quadratic_form,'pe')
     eff=X(:,1:N)*coeff(1:N);
     if K>0
        A_b=[X(:,1:N)'*(eff-mean(eff)); sparse(J-1,1); zeros(K,1)];
     end
     if K==0
        A_b=[X(:,1:N)'*(eff-mean(eff)); sparse(J-1,1)];
     end
               
  end
               
  if strcmp(type_quadratic_form,'cov')
     pe=X(:,1:N)*coeff(1:N,1);
     fe=X(:,N+1:N+J-1)*coeff(N+1:N+J-1,1);
     
     if K>0
        A_b=[0.5*X(:,1:N)'*(fe-mean(fe)); 0.5*X(:,N+1:N+J-1)'*(pe-mean(pe)); zeros(K,1)];
     end
     
     if K==0
        A_b=[0.5*X(:,1:N)'*(fe-mean(fe)); 0.5*X(:,N+1:N+J-1)'*(pe-mean(pe))];
     end   
  end
  
%This completes the right part
 C_vec = construc_W(veC,X,xx,Lchol,A_b,Lambda_B,I_Lambda_P,L_P,res_h);

end

