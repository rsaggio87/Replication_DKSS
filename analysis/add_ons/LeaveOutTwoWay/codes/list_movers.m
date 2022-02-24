function list = list_movers(firmid,id)

%find movers
    J=max(firmid);
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
    
%id and firmids associated with movers only
    id_movers=id(move);
    firmid_mover=firmid(move);

%need to normalize id_movers and keep a dictionary. 
    id_movers_orig=id_movers;
    [~,~,id_movers]=unique(id_movers);
   
    
%unique pairs    
    list=[id_movers firmid_mover id_movers_orig];
    list=unique(list,'rows');
    
 end   