% this script describes some of the interactions 
% that mediate the effects of Abeta on synapses


% the first three blocks should be commented out
% before running absynCMB

% run absynINI to initialize the simulation
absynINI

% % the inputs are presyanptic activity and Abeta
% % preSYN is an element of [0 1 2 3], Abeta of [0 1]
% preSYN = 0;
% Abeta = 0;
% % 
% % set combination of compounds
% % for receptors
% AChRnorm = 0;
% mGRblock = 0;
% TrkBnorm = 0;
% % for kinases and associated
% ACact = 0;
% GSK3block = 0;
% PDEblock = 0;
% PKCact = 0;
% % for phosphatases and associated
% PP1block = 0;
% PP2Bblock = 0;
% proACT = 0;
% % block for PKC is experimental only
% PKCblock = 0;

% set constants
ACthr = 7;
BDNFthr = 2;
CaMKthr= 9;
CaSAT = 20;
PKCthr = 8;
PP2Bthr = 5;

% set loop parameters
countLIM = 20; % enter count limit
counter = 0; % zero counter
upFlag = 1; % set upFlag to one
CaRec = zeros(1, countLIM); % set a Ca record

% update transmitters and mediators
Glu = preSYN;
ACh = preSYN;
if preSYN > BDNFthr
    BDNF = 1; 
else
    BDNF = 0;
end

% enter while loop
while upFlag ~= 0
    counter = counter + 1;
    if counter >= countLIM, break, end
    upFlag = 0;
    
    if Glu==0
        NMDARhld = 0;
    elseif Glu > 0
        NMDARhld = NR1 + NR2 + Glu;
    end
    if NMDAR ~= NMDARhld
        NMDAR = NMDARhld;
        upFlag = 1;
    end
    
    if Glu==0 | mGRblock==1
        mGluR5hld = 0;
    elseif Glu > 0 && mGRblock==0
        mGluR5hld = 1;
    end
    if mGluR5 ~= mGluR5hld
        mGluR5 = mGluR5hld;
        upFlag = 1;
    end
    
    if Abeta==0 | AChRnorm==1
        nAChRhld = ACh;
    elseif Abeta==1 && AChRnorm==0 && ACh==0   
        nAChRhld = 1;
    elseif Abeta==1 && AChRnorm==0 && ACh > 0;
        nAChRhld = 2;
    end
    if nAChR ~= nAChRhld
        nAChR = nAChRhld;
        upFlag = 1;
    end
    
    if TrkB ~= BDNF
        TrkB = BDNF;
        upFlag = 1;
    end
    
    if Abeta==0 | TrkBnorm==1
        Shchld = TrkB;
    elseif Abeta==1 && TrkBnorm==0
        Shchld = 0;
    end
    if Shc ~= Shchld
        Shc = Shchld;
        upFlag = 1;
    end
    
    if SOS ~= Shc
        SOS = Shc;
        upFlag = 1;
    end
    
    if Abeta==0 | TrkBnorm==1
        IRShld = TrkB;
    elseif Abeta==1 && TrkBnorm==0
        IRShld = 0;
    end
    if IRS ~= IRShld
        IRS = IRShld;
        upFlag = 1;
    end
    
    if Ras ~= SOS
        Ras = SOS;
        upFlag = 1;
    end
    
    if Ras==1 | IRS==1
        PI3Khld = 1;
    else
        PI3Khld = 0;
    end
    if PI3K ~= PI3Khld
        PI3K = PI3Khld;
        upFlag = 1;
    end
    
    if PDK1 ~= PI3K
        PDK1 = PI3K;
        upFlag = 1;
    end
    
    if Akt ~= PDK1
        Akt = PDK1;
        upFlag = 1;
    end
    
    if Gq ~= mGluR5
        Gq = mGluR5;
        upFlag = 1;
    end
    
    if PLC ~= Gq
        PLC = Gq;
        upFlag = 1;
    end
    
    if PIP2 > 0 && IP3 ~= PLC
        IP3 = PLC;
        upFlag = 1;
    end
    
    if PIP2 > 0 && DAG ~= PLC
        DAG = PLC;
        upFlag = 1;
    end
    
    if IP3R1 ~= IP3
        IP3R1 = IP3;
        upFlag = 1;
    end
    
    Cahld = NMDAR + nAChR + IP3R1;
    if Cahld > CaSAT, Cahld = CaSAT; end
    if Ca ~= Cahld
        Ca = Cahld;
        upFlag = 1;
        CaRec(counter) = Ca;
    end

    if Ca <= PKCthr && PKCact==0
        PKChld = 0;
    elseif Ca > PKCthr && PKCact==0 && PKCblock==0  
        PKChld = 1;
    elseif PKCblock==1
        PKChld = 0;
    elseif PKCblock==0 && PKCact==1
        PKChld = 1;
    end
    if PKC ~= PKChld
        PKC = PKChld;
        upFlag = 1;
    end
    
    if CaM ~= Ca
        CaM = Ca;
        upFlag = 1;
    end
    
    if CaM <= ACthr && ACact==0
        AChld = 0;
    elseif CaM > ACthr | ACact==1
        AChld = 1; 
    end
    if AC~= AChld
        AC = AChld;
        upFlag = 1;
    end
    
    if PDE ~= 1 - PDEblock
        PDE = 1 - PDEblock; 
        upFlag = 1;
    end
    
    if AC == 0
        cAMPhld = 0;
    elseif AC > 0
        cAMPhld = 2 - PDE + AC;
    end
    if cAMP ~= cAMPhld
        cAMP = cAMPhld;
        upFlag = 1;
    end
    
    if PKA ~= cAMP
        PKA = cAMP;
        upFlag = 1;
    end
    
    if PP2Bblock==0 && CaM > PP2Bthr
        PP2Bhld = 1; 
    elseif CaM <= PP2Bthr | PP2Bblock==1
        PP2Bhld = 0;
    end
    if PP2B ~= PP2Bhld
        PP2B = PP2Bhld;
        upFlag = 1;
    end
    
    if PKA == 0 && PP2B == 0
        PPI1hld = 1;
    elseif PKA == 0 && PP2B == 1
        PPI1hld = 0;
    elseif PKA > 0 
        PPI1hld = 1;
    end
    if PPI1 ~= PPI1hld
        PPI1 = PPI1hld;
        upFlag = 1;
    end
    
    if Abeta==0 && proACT==0
        prothld = 1;
    elseif Abeta==1 && proACT==0
        prothld = 0;
    elseif proACT==1
        prothld = 2;
    end
    if proteosome ~= prothld
        proteosome = prothld;
        upFlag = 1;
    end
        
    if proteosome == 2
        STEPhld = 0;
    elseif proteosome == 1
        STEPhld = max(PP2B - PKA, 0);
    elseif proteosome == 0 && PP2B > 0
        STEPhld = max(PP2B + 1 - PKA, 0);
    end
    if STEP ~= STEPhld
        STEP = STEPhld;
        upFlag = 1;
    end
    
    if Pyk2 ~= max(PKC - STEP, 0)
        Pyk2 = max(PKC - STEP, 0);
        upFlag = 1;
    end
    
    if Src ~= Pyk2
        Src = Pyk2;
        upFlag = 1;
    end
    
    if Src == 0
        NR2hld = 1;
    elseif Src > 0,
        NR2hld = 2;
    end
    if NR2 ~= NR2hld
        NR2 = NR2hld;
        upFlag = 1;
    end
    
    if PP1 > 0 && CaM <= CaMKthr
        CaMKIIhold = 0;
    elseif CaM > CaMKthr
        CaMKIIhld = 1; 
    end
    if CaMKII ~= CaMKIIhld
        CaMKII = CaMKIIhld;
        upFlag = 1;
    end
    
    if PP1block==0 && PPI1==0 && CaMKII==0 && I2==0 
        PP1hld = 1;
    else
        PP1hld = 0;
    end
    if PP1 ~= PP1hld
        PP1 = PP1hld;
        upFlag = 1;
    end

    if GSK3block==1
        GSK3hld = 0;
    elseif Akt==0 && PP1==0 
        GSK3hld = 1;
    elseif Akt==0 && PP1==1 
        GSK3hld = 2;
    elseif Akt==1 && PKA==2 && PKC==1 && CaMKII==1
        GSK3hld = 0;
    elseif Akt==1
        GSK3hld = max(1 + PP1 - Akt - min(PKA, 1), 0);
    end
    if GSK3 ~= GSK3hld
        GSK3 = GSK3hld;
        upFlag = 1;
    end    
           
    if GSK3 > 0
        I2hld = 0;
    else
        I2hld = 1;
    end
    if I2 ~= I2hld
        I2 = I2hld;
        upFlag = 1;
    end

    AR1hld = max(1 + (CaMKII + PKC + max(PKA-1, 0))-(PP1 + PP2B), 0);
    if AR1 ~= AR1hld
        AR1 = AR1hld;
        upFlag = 1;
    end

    if AR2 ~= max(1 - STEP, 0)
        AR2 = max(1 - STEP, 0);
        upFlag = 1;
    end
    
    if KLC2 ~= max(2 - GSK3, 0)
        KLC2 = max(2 - GSK3, 0);
        upFlag = 1;
    end
    
    if AMPAR ~= max(AR1 + AR2 + KLC2, 0)
        AMPAR = max(AR1 + AR2 + KLC2, 0);
        upFlag = 1;
    end
    
end % end while loop


return

% display values
sprintf('\n \n \n *** NEW RESULTS ***')
Ca
% AChRnorm
% mGRblock
% TrkBnorm
% PDEblock
% PP2B
% PP1
% GSK3
% PKA
% PKC
% CaMKII
% NMDAR
% I2
% proteosome
% STEP
AMPAR
counter


