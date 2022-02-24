				data=importdata('../tmp/RovigoCS_Leave2.csv');
				id=data(:,4);
				firmid=data(:,2);
				y=data(:,3);
				year=data(:,end);
				controls=year==2001;
				clear data
				
				[~,~,id]=unique(id);
				[~,~,firmid]=unique(firmid);
				
				gcs = [NaN; id(1:end-1)];
				gcs = id~=gcs;
				lagfirmid=[NaN; firmid(1:end-1)];
				lagfirmid(gcs==1)=NaN; %%first obs for each worker
				stayer=(firmid==lagfirmid);
				stayer(gcs==1)=1;
				stayer=accumarray(id,stayer);
				T=accumarray(id,1);
				stayer=T==stayer;
				movers=stayer~=1;
				movers=movers(id);
				T=T(id);
				id_movers=id(movers);
				[~,~,n]=unique(id_movers);
				Nmovers=max(n);
				NT=size(y,1);
				D=sparse(1:NT,id',1);
				N=size(D,2);
				F=sparse(1:NT,firmid',1);
				J=size(F,2);
				stayer=stayer(id);
				
				X=[D -F];
				xx=X'*X;
				for ii=1:100
				xtilde= xx\X(ii,:)';
				aux_right=xtilde(N+1:end);
                aux_left=xtilde(1:N);
				%COV=cov(D*aux_left,F*aux_right);
				COV=cov(F*aux_right,F*aux_right);
				num2str([COV(1,2) stayer(ii)])
				end
				