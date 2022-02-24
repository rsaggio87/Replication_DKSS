function TABELLA=table_2(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename)

%Run KSS
[one two three] = leave_out_COMPLETE(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename);
						 
%Load the auxiliary file for table
s=[filename '_completed'];
load(s)
pe=D*ahat;

%Export the table
TABELLA(1) = var(y);

TABELLA(2) = sigma_2_psi_AKM;
TABELLA(3) = var_corrected_fe;
TABELLA(4) = sigma2_psi;

TABELLA(5) = sigma_alpha_psi_AKM;
TABELLA(6) = var_corrected_cov;
TABELLA(7) = sigma_psi_alpha;

TABELLA(8) = sigma_2_alpha_AKM;
TABELLA(9) = var_corrected_pe;
TABELLA(10) = sigma2_alpha;

TABELLA(11) = corr(fe,pe);
TABELLA(12) = var_corrected_cov/(sqrt(var_corrected_fe)*sqrt(var_corrected_pe));
TABELLA(13) = sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha));

TABELLA(14) = R2;
TABELLA(15) = adjR2;
TABELLA(16) = explained_var_leave_out;


%Export to .csv to visualize network
if 0 == 1
			gcs = [NaN; id(1:end-1)];
			gcs = id~=gcs;
			lagfirmid=[NaN; firmid_old(1:end-1)];
			lagfirmid(gcs==1)=NaN; %%first obs for each worker
			sel=~isnan(lagfirmid);
			out=[id_old(sel) lagfirmid(sel) firmid_old(sel) controls(sel) x1bar_all(sel,1:2)]; %%all the moves observed from one period to the next
			s=[filename 'MIKKEL.csv'];
			dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
end


end


