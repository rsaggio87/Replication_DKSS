function sigma_predict= llr_fit(Pii,Bii,sigma_i,subsample_llr_fit,K,movers,T,KGrid,aa)

if nargin <= 7
    KGrid = 1000;
    aa = 1; 
end

if nargin <= 8
    aa = 1; 
end

if size(KGrid,1) == 0
   KGrid = 1000;     
end

%% Set up dimensions
NT=size(sigma_i,1);
maxT=max(T);

Nmovers=sum(movers);
%% Run LLR
hBest=min(aa/NT^(1/3),1);
if subsample_llr_fit == 0
    f = fit([Pii Bii],sigma_i,'lowess','Normalize','on','Span',hBest);
    sigma_predict = feval(f,[Pii, Bii]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 1 && K==0
    sigma_predict=zeros(NT,1);
    for ti=1:maxT
            sel=(T==ti & ~movers);
            sigma_predict(sel)=mean(sigma_i(sel));
    end
    if Nmovers > 0
    sigma_use=sigma_i(movers);
    Xuse=[Bii(movers) Pii(movers)];
    cellS=size(Xuse,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest);
    sigma_predict(movers) = feval(f,[Bii(movers), Pii(movers)]);
    end
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 1 && K>0
    sigma_predict=zeros(NT,1);
    %LLR on Stayers.
    sigma_use=sigma_i(~movers);
    Xuse=[Bii(~movers) Pii(~movers)];
    cellS=size(Xuse,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest);
    sigma_predict(~movers) = feval(f,[Bii(~movers), Pii(~movers)]);
    %LLR on Movers.
    if Nmovers>0
    Xuse=[Bii(~movers) Pii(~movers)];
    cellS=size(Xuse,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest);
    sigma_predict(movers) = feval(f,[Bii(movers), Pii(movers)]);
    end
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
    
end

if subsample_llr_fit == 2 && K==0
    sigma_stayers=zeros(NT,1);
    for ti=1:maxT
            sel=(T==ti & ~movers);
            sigma_stayers=sigma_stayers+sparse(find(sel),1,mean(sigma_i(sel)),NT,1);
    end
    sigma_movers=0;
    if Nmovers>0
    Bii=Bii(movers);
    Pii=Pii(movers);
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_i(movers),[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_movers=sparse(find(movers),1,feval(f,[Bii, Pii]),NT,1);
    end
    sigma_predict = sigma_movers + sigma_stayers;
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 2 && K>0
    Bii_old=Bii;
    Pii_old=Pii;
    %Binned LLR on stayers.
    Bii=Bii_old(~movers);
    Pii=Pii_old(~movers);
    sigma_use=sigma_i(~movers);
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_use,[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_stayers=sparse(find(~movers),1,feval(f,[Bii, Pii]),NT,1);
    %Binned LLR on movers.
    sigma_movers=0;
    if Nmovers > 0
    Bii=Bii_old(movers);
    Pii=Pii_old(movers);
    sigma_use=sigma_i(movers);
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_use,[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_movers=sparse(find(movers),1,feval(f,[Bii, Pii]),NT,1);
    end
    sigma_predict = sigma_stayers + sigma_movers;
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 3
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_i,[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=min(aa/cellS^(1/3),1);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_predict= feval(f,[Bii, Pii]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 4
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    sigma_use=accumarray(qq,sigma_i,[],@mean);
    sigma_predict=sigma_use(qq);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

end    


