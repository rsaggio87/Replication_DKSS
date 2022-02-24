function TABELLA=table_1(y,id,firmid,controls,city)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes as input the person year observations from a particular
%region and then performs the following tasks:

%IMPORTANT: original src data must be sorted by id year (xtset id year)

% 1. Finds and summarizes the largest connected set using routine in CHK (2013).
% 2. Proceeds to find leave out largest connected set, see Appendix B.
% 3. Proceed to find the largest leave 2 out connected set, see Appendix B.

%Along the way, the code will output important information about each sample.
%Information on the maximum leverage will be computed while performing the analysis on
%the variance decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%path to matlab bgl
TABELLA=zeros(15,1);

%% STEP 1: FIND LARGEST CONNECTED SET
%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker

%Largest connected set.
[y,id,firmid,id_old,firmid_old,controls] = connected_set(y,id,firmid,lagfirmid,controls);

%Find movers
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
Nmovers=sum(movers);
J=max(firmid);

clear gcs stayer movers T lagfirmid

s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Results on Largest Connected Set '];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(max(Nmovers))];
disp(s);
s=['# of Firms: ' num2str(J)];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);


%SAVE LARGEST CONNECTED SET
TABELLA(1)=size(y,1);
TABELLA(2)=max(Nmovers);
TABELLA(3)=J;
TABELLA(4)=mean(y);
TABELLA(5)=var(y);

%% STEP 2: LEAVE ONE OUT CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave one out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic
[y,firmid,id,id_old,firmid_old,controls]= pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);
disp('Time to find leave one out largest connected set')
toc
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%%%Drop stayers with a single person year observation
T=accumarray(id,1);
T=T(id);
sel=T>1;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);

%Resetting ids one last time.
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n; 

%Find movers
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
Nmovers=sum(movers);
movers=movers(id);
J=max(firmid);
clear gcs T lagfirmid


%% STEP 3: EXTRA PRUNING FOR ROVIGO ONLY
%We provide one extra pruning for Rovigo to make sure that conditions of
%Theorem 1 are satisfied. In particular, we identify a mover with the
%highest Lindeberg condition and remove the stayers from the two frms this 
%individual moved between. 

%To understand how we identified this particular influential observation,
%comment the following lines and proceed to calculations of the Lindeberg
%condition on the unpruned leave out sample of Rovigo. This will show that 
%mover associated with the given py observation is the one with the
%highest Lindeberg condition. 

if strcmp(city,'Rovigo')
firms_bad=[129;151];
sel=ismember(firmid,firms_bad);
sel=and(sel,~movers);
sel=~sel; 
y=y(sel);
id=id(sel);
firmid=firmid(sel);
controls=controls(sel,:);  
id_old=id_old(sel);
firmid_old=firmid_old(sel);
end

%% STEP 4: SUMMARIZE LEAVE OUT CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(J)];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%SAVE LEAVE OUT LARGEST CONNECTED SET
TABELLA(6)=size(y,1);
TABELLA(7)=max(Nmovers);
TABELLA(8)=J;
TABELLA(9)=mean(y);
TABELLA(10)=var(y);
out=[id,firmid,y,id_old,firmid_old,controls];
s=['tmp/' city 'CS_Leave1.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 



%% STEP 5: LEAVE TWO OUT CONNECTED SET
tic
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave two out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic
[id,firmid,y,controls,id_old,firmid_old] = pruning_leave2out(id,firmid,y,controls,id_old,firmid_old);
disp('Time to find leave one out largest connected set')
toc
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%%%Drop stayers with a single person year observation
T=accumarray(id,1);
T=T(id);
sel=T>1;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);

%Resetting ids one last time.
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n; 

%Find movers
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
Nmovers=sum(movers);
movers=movers(id);
J=max(firmid);
clear gcs T lagfirmid

%% STEP 6: SUMMARIZE LEAVE TWO OUT CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['Info on the leave two out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(J)];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%SAVE LEAVE OUT LARGEST CONNECTED SET
TABELLA(11)=size(y,1);
TABELLA(12)=max(Nmovers);
TABELLA(13)=J;
TABELLA(14)=mean(y);
TABELLA(15)=var(y);
out=[id,firmid,y,id_old,firmid_old,controls];
s=['tmp/' city 'CS_Leave2.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 


end


