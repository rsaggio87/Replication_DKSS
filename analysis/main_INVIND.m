%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This m-file computes the DWL model of Di Addario, Kline, Saggio and
%Soelvsten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %External Packages                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
path(path,'add_ons/matlab_bgl/'); %note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
path(path,'add_ons/LeaveOutTwoWay/codes'); %this contains the main LeaveOut Routines.
path(path,'codes');
warning('off', 'all') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Table 1,2,3,4(column 1) 5 and 6.                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp = 1:3
%Read options
			type_of_algorithm = 'JLL'
			scale			  = 250;
			random_sample	  = 100;
			pool			  = pp;
            J_J_STRINGA       = '_INVIND_data_1990_2016_2005_2016_RS_build'

%Raw data
			namesrc	=['../build/src/RAKIM/INVIND_RAKIM' J_J_STRINGA '.csv']
			
%Log File

			if pool == 1
			
			poolSTRINGA=['POOLED_SAMPLE_' J_J_STRINGA]
			 
			end
			
			if pool == 2
			
			poolSTRINGA=['MALE_SAMPLE_' J_J_STRINGA]
			
			end
			
			if pool == 3
			
			poolSTRINGA=['FEMALE_SAMPLE_' J_J_STRINGA]
			
            end
			placelog='logs/';
            namelog=['BUILD_INVIND' poolSTRINGA '_RAKIM_' type_of_algorithm '_' num2str(random_sample)];
			logname=['logs/' namelog  '.log'];
			system(['rm ' logname])
			diary(logname)

			
%Import data
			data=importdata(namesrc);
			
%Read specification			
			
			if pool == 1
			sel 	=data(:,1)>0;
			end
			
			if pool == 2
			sel 	=data(:,7)==0; %male only
			end
			
			if pool == 3
			sel 	=data(:,7)==1; %female only
			end
			
			data				= data(sel,:);
			
%Read variables			
			id					= data(:,1);
			firmid				= data(:,2);
			lagfirmid			= data(:,3);
			year				= data(:,4);
			y					= data(:,5);
			age					= data(:,6);
			female				= data(:,7);
			[~,~,controls]  	= unique(year);
			controls 	    	= [sparse((1:size(y,1))',controls',1,size(y,1),max(controls))];
			controls(:,end) 	= [];
			controls			= [controls ((age-40)/40).^2 ((age-40)/40).^3];
						
%Renormalize worker ids		
			[~,~,id]			= unique(id);
					
													
%Random sample?
if random_sample < 100
			rng(1234)		
			idsel				= randi(max(id),round((random_sample/100)*max(id)),1);
			sel					= ismember(id,idsel);
			id					= id(sel);
			firmid				= firmid(sel);
			lagfirmid			= lagfirmid(sel);
			y					= y(sel);
			controls			= controls(sel,:);
			female				= female(sel);
end

%STEP 1: FIND IDENTIFIED SET (THIS TAKES A LONG-TIME)
			identified			= kss_rakim_build(id,firmid,lagfirmid,y,controls,female,scale,poolSTRINGA);
	
%STEP 2: ESTIMATE THE MODEL (log-file prints Table 1, 2, 3, 5 and 6)
			s					= ['../build/src/RAKIM/PRUNED_' poolSTRINGA '.csv']
			data				= importdata(s);
			y					= data(:,1);
			id					= data(:,2);
			firmid_orig			= data(:,3);
			lagfirmid_orig		= data(:,4);
			firmid				= data(:,5);
			lagfirmid			= data(:,6);
			female				= data(:,7);
			controls			= data(:,8:end);
			clear data					
			theta				= kss_rakim_levels(y,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,female,scale,poolSTRINGA);								
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %Table 4               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cc=2:4 %looping over columns 2,3 and 4 of Table 4
    		
            if 	cc == 2
			J_J_STRINGA = '_DOMINANT_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build'
			end	
			
			if 	cc == 3
			J_J_STRINGA = '_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build'
			end				

			if 	cc == 4
			J_J_STRINGA = '_WITH_STAYERS_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build'
            end	
            
            poolSTRINGA=['POOLED_SAMPLE_' J_J_STRINGA]
            
 %Read
            placelog='logs/';
            namelog=['BUILD_INVIND' poolSTRINGA '_RAKIM_'];
			logname=['logs/' namelog  '.log'];
			system(['rm ' logname])
			diary(logname)

			
%Import data
            namesrc	=['../build/src/RAKIM/INVIND_RAKIM' J_J_STRINGA '.csv']
			data=importdata(namesrc);
			
%Read specification			
			sel                     = data(:,1)>0;
			data                    = data(sel,:);
			
%Read variables			
			id                      = data(:,1);
			firmid                  = data(:,2);
			lagfirmid               = data(:,3);
			year                    = data(:,4);
			y                       = data(:,5);
			age                     = data(:,6);
			female                  = data(:,7);
			[~,~,controls]          = unique(year);
			controls                = [sparse((1:size(y,1))',controls',1,size(y,1),max(controls))];
			controls(:,end)         = [];
			controls                = [controls ((age-40)/40).^2 ((age-40)/40).^3];
						
%Renormalize worker ids		
			[~,~,id]                = unique(id);
            
%Options to run KSS  
            resid_controls			= 1;
			andrews_estimates		= 0;
			eigen_diagno			= 0;
			subsample_llr_fit		= 0;
			restrict_movers 		= 0;
			type_of_algorithm		='JLL'
			epsilon					= 0.01;
			filename				= ['/scratch/public/leave_out/KSS_' poolSTRINGA];
			do_SE					= 0;
			cd 'add_ons/LeaveOutTwoWay/codes'   
			MakeCMG;
			cd '../../'	
%KSS					
			leave_out_level 		= 'matches'
			[sigma2_psi] 	    	= leave_out_COMPLETE(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename);
            
end
diary off			 
delete(gcp('nocreate'))
