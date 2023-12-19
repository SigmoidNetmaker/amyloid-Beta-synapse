% this script finds specific response sets 
% from repARRAY as computed using absynCMB
% which runs absyn

% the structure of the arrays in absynCMB must
% be taken into accound for all of the follwing

% enter an ordered response set to be found
% probe = [3 3 0 6];
% then paste the following at beginning of loop
%     if repARRAY(i,1+nDRUGS+1:1+nDRUGS+4)==probe

% or devise some other set of conditionals instead
% and paste those at the beginning of the loop
% various normalizations of plasticity
%     if repARRAY(i,1+nDRUGS+1) == repARRAY(i,1+nDRUGS+2) && ...
%        repARRAY(i,1+nDRUGS+1)  > repARRAY(i,1+nDRUGS+3) && ...
%        repARRAY(i,1+nDRUGS+1)  < repARRAY(i,1+nDRUGS+4) && ...
%        repARRAY(i,2) == 0

%     if repARRAY(i,1+nDRUGS+1) == repARRAY(i,1+nDRUGS+2) && ...
%        repARRAY(i,1+nDRUGS+1)  > repARRAY(i,1+nDRUGS+3)

%     if repARRAY(i,1+nDRUGS+1) == repARRAY(i,1+nDRUGS+2) && ...
%        repARRAY(i,1+nDRUGS+1)  > repARRAY(i,1+nDRUGS+3)

%     if repARRAY(i,1+nDRUGS+2) > repARRAY(i,1+nDRUGS+3) && ...
%        repARRAY(i,1+nDRUGS+2) < repARRAY(i,1+nDRUGS+4) && ...
%        repARRAY(i,2) ~= 1

%     if repARRAY(i,1+nDRUGS+1) < repARRAY(i,1+nDRUGS+4) && ...
%        repARRAY(i,2)==0 && repARRAY(i,4)==0 && repARRAY(i,10)==0

%     if repARRAY(i,1+nDRUGS+1) < repARRAY(i,1+nDRUGS+4) && ...
%        repARRAY(i,2)==0 && repARRAY(i,3)==0 && repARRAY(i,4)==0

% improvement in LTD without receptor compounds
%     if repARRAY(i,1+nDRUGS+2) > repARRAY(i,1+nDRUGS+3) && ...
%             repARRAY(i,2)==0 && repARRAY(i,4)==0 && repARRAY(i,3)==0 

% improvement on AChRnorm to prefect LTD and LTP with non-receptor
%     if repARRAY(i,1+nDRUGS+1)==3 && repARRAY(i,1+nDRUGS+2)==3 && ...
%         repARRAY(i,1+nDRUGS+3)==0 && repARRAY(i,1+nDRUGS+4) > 5 && ...
%         repARRAY(i,2)==1 && repARRAY(i,3)==0 && repARRAY(i,4)==0

% improvement on mGRblock to prefect LTD and LTP with non-receptor
%     if repARRAY(i,1+nDRUGS+1)==3 && repARRAY(i,1+nDRUGS+2)==3 && ...
%         repARRAY(i,1+nDRUGS+3)==0 && repARRAY(i,1+nDRUGS+4) > 3 && ...
%         repARRAY(i,2)==0 && repARRAY(i,3)==1 && repARRAY(i,4)==0

% any LTD and LTP improvement on mGRblock with non-receptor
%     if repARRAY(i,1+nDRUGS+2) > repARRAY(i,1+nDRUGS+3) && ...
%         repARRAY(i,1+nDRUGS+4) > repARRAY(i,1+nDRUGS+1) && ...
%         repARRAY(i,2)==0 && repARRAY(i,3)==1 && repARRAY(i,4)==0

% any LTD and LTP improvement on TrkBnorm with non-receptor
%     if repARRAY(i,1+nDRUGS+2) > repARRAY(i,1+nDRUGS+3) && ...
%         repARRAY(i,1+nDRUGS+4) > repARRAY(i,1+nDRUGS+1) && ...
%         repARRAY(i,2)==0 && repARRAY(i,3)==0 && repARRAY(i,4)==1

% any LTD and LTP improvement on generic non-receptor with non-receptor
%     if repARRAY(i,1+nDRUGS+2) > repARRAY(i,1+nDRUGS+3) && ...
%         repARRAY(i,1+nDRUGS+4) > repARRAY(i,1+nDRUGS+1) && ...
%         repARRAY(i,2)==0 && repARRAY(i,3)==0 && repARRAY(i,4)==0 && ...
%         repARRAY(i,?)==1 % where ? is the non-receptor in question

% any LTP improvement on generic non-receptor with non-receptor
%     if repARRAY(i,1+nDRUGS+4) > repARRAY(i,1+nDRUGS+1) && ...
%         repARRAY(i,2)==0 && repARRAY(i,3)==0 && repARRAY(i,4)==0 && ...
%         repARRAY(i,?)==1 % where ? is the non-receptor in question

% generate the truth vector
truthVEC = zeros(nCMBS,1);
for i = 1:nCMBS
    if repARRAY(i,1+nDRUGS+4) > repARRAY(i,1+nDRUGS+1) ...
        && repARRAY(i,2)==0 && repARRAY(i,3)==0 && repARRAY(i,4)==0 ...
        && repARRAY(i,11)==1
        truthVEC(i) = 1;
    end
end

% find combos matching probe
hits=find(truthVEC);
hitREP=repARRAY(hits,:)
% hitRES=resARRAY(hits,:)
% hitCMBS=cmbARRAY(hits,:)




    