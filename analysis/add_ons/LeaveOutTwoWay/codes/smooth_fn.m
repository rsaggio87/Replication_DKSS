function sigma_predict= smooth_fn(elist,Lambda_P,Lambda_B,year,y,eta_h,movers,K,T,subsample_llr_fit,leave_out_level,restrict_movers,type_quadratic_form,bw)
%This function performs the smoothing step required to conduct inference 
%on variance components. This function can account for the presence of serial 
%correlation within match, as indicated by the user when specifying 
%leave_out_level='match'.

    tic
    NT=size(Lambda_P,1);

    rows=elist(:,1);
    columns=elist(:,2);
    linearInd = sub2ind(size(Lambda_P),rows,columns);
    NSTAR=size(rows,1);

    Pii = Lambda_P(linearInd);
    Bii = Lambda_B(linearInd);
    sigma_i = 0.5*(y(rows).*eta_h(columns)+y(columns).*eta_h(rows));

    movers_cl = movers(rows);
    T_cl = T(rows);

    if  strcmp(leave_out_level,'obs')
        sigma_predict = llr_fit(Pii,Bii,sigma_i,subsample_llr_fit,K,movers_cl,T_cl,[],bw);      
    end

    if  ~strcmp(leave_out_level,'obs') 
                    year_r = year(rows);
                    year_c = year(columns);
                    [~,~,cmb]=unique([year_r year_c],'rows');
                    Ncombo=max(cmb);
                    sigma_predict=zeros(NSTAR,1);
                    parfor s = 1:Ncombo
                           sel = (cmb == s);
                           Pii_use = Pii(sel);
                           Bii_use = Bii(sel);
                           sig_use = sigma_i(sel);
                           mov_use = movers_cl(sel);
                           Tto_use = T_cl(sel);
                           CONTO = sum(mov_use);

                           if CONTO >= 1000
                               sigma_aux=llr_fit(Pii_use,Bii_use,sig_use,subsample_llr_fit,K,mov_use,Tto_use,[],bw);    
                           end


                           %all the steps below account for situations where
                           %the user wants to run a binned up version of the
                           %smoothing step but it has too few data points.

                           if CONTO <= 1000 && subsample_llr_fit == 1
                                sigma_aux=llr_fit(Pii_use,Bii_use,sig_use,1,K,mov_use,Tto_use,[],bw);  %no action required.  
                           end

                           if CONTO <= 1000 && subsample_llr_fit == 2
                                sigma_aux=llr_fit(Pii_use,Bii_use,sig_use,1,K,mov_use,Tto_use,[],bw);    
                           end
                           
                           if CONTO <= 1000 && subsample_llr_fit == 3
                                sigma_aux=llr_fit(Pii_use,Bii_use,sig_use,0,K,mov_use,Tto_use,[],bw);    
                           end
                           
                           if CONTO <= 1000 && subsample_llr_fit == 4 && CONTO>10
                                sigma_aux=llr_fit(Pii_use,Bii_use,sig_use,4,K,mov_use,Tto_use,10,bw);    
                           end
                           
                           if CONTO <=100 && subsample_llr_fit == 4
                           error('Too few data points in bucket where py observations correspond to')
                           [mean(year_r(sel)) mean(year_c(sel))]
                           end
                           %collect
                           sigma_predict=sigma_predict+sparse(find(sel),1,sigma_aux,NSTAR,1);
                   end
    end

    %Bring back to NT X NT matrix.
    if strcmp(leave_out_level,'matches') && restrict_movers == 0 && K == 0 && ~strcmp(type_quadratic_form,'pe')
        rows=rows(movers_cl);
        columns=columns(movers_cl);
        sigma_predict=sigma_predict(movers_cl);
    %I do this because I know that for cov(pe,fe) and var(fe) the matrix C has
    %all zeros for stayers when there are no controls in the model.
    end  
    
    
    sigma_predict = sparse(rows,columns,sigma_predict,NT,NT);
    sigma_predict = sigma_predict + triu(sigma_predict,1)'; %make it symmetric.
    disp('Time to Perform computation of the standard error')
    toc
end

