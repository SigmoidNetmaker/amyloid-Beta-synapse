# amyloid-Beta-synapse
MATLAB and Maude code for representing, simulating, and analyzing the molecular interactions by which amyloid-Beta affects synapses, specifically long-term potentiation (LTP). As a declarative language, Maude was used to specify the cellular and molecular interactions that mediate synaptic plasticity, and then analyze that specification to yeild insight into amyloid-Beta pathology. MATLAB was used both to cross-check the Maude specification and for more computationally intensive procedures.  

These are the Maude files used to specify and analyze synaptic amyloid-Beta interactions:  

absyn.txt -- this module specifies the cell-molecular interactions by which amyloid-Beta affects LTP  
mc-absyn.txt -- this module sets up model-checking of the amyloid-Beta-LTP specification  
absyn-check.txt -- this command model-checks the amyloid-Beta-LTP specification  
absyn-rew.txt -- this command executes rewrites of the amyloid-Beta-LTP specification  
absyn-search.txt -- this command executes a state-space search of the amyloid-Beta-LTP specification  

These are the MATLAB files used to cross-check Maude and for other procedures:  

absyn.m -- this script describes the same interactions as the Maude module absyn.txt shown above  


