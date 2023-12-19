% this script explores drug combinations in absyn.m
% the first three code blocks in absyn.m need to be
% out commented before this script is run

% enter number of drugs 
nDRUGS = 10;
nCMBS = 2^nDRUGS;

% generate combo number vector and 
% binary array for combinations
cmbNUM = (1:2^nDRUGS)';
cmbARRAY = rem(floor((cmbNUM - 1)*pow2(-(nDRUGS-1):0)),2);

% set up report array and results array
repARRAY = [cmbNUM cmbARRAY zeros(nCMBS,4)];
resARRAY = [cmbNUM zeros(nCMBS,4)];

% evaluate results of combinations with Abeta=1
% note that PKCblock also is zero
Abeta = 1;
PKCblock = 0;
for cmb = 1:nCMBS
    for act = 1:4
        AChRnorm = cmbARRAY(cmb,1);
        mGRblock = cmbARRAY(cmb,2);
        TrkBnorm = cmbARRAY(cmb,3);
        ACact = cmbARRAY(cmb,4);
        GSK3block = cmbARRAY(cmb,5);
        PDEblock = cmbARRAY(cmb,6);
        PKCact = cmbARRAY(cmb,7);
        PP1block = cmbARRAY(cmb,8);
        PP2Bblock = cmbARRAY(cmb,9);
        proACT = cmbARRAY(cmb,10);
        preSYN = act - 1;
        absynINI
        absyn
        repARRAY(cmb,1+nDRUGS+act) = AMPAR;
        resARRAY(cmb,1+act) = AMPAR;
    end
end









