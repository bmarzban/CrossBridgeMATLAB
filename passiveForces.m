function [PF,dPF] = passiveForces(SL,SL_rest,kpe1)
% Passive Force formulation adopted from Rice etal Rice etal (Biophys J. 2008 Sep;95(5):2368-90)
% Function to calculate the passive forces of the muscle
% The passive force is a function of the sarcomere length
beta = 3.5; % um
PCon_titin = kpe1*beta*0.002; %
PExp_titin = 10; %(um-1)
SL_collagen = 2.25; %(um)
PCon_collagen = kpe1*beta*0.02; %
PExp_collagen  = 70; %(um-1)

% Passive forces: Trabeculae 
PF_titin = sign(SL-SL_rest)*PCon_titin*exp(abs(SL-SL_rest)*PExp_titin-1);
dPF_titin = (PCon_titin*PExp_titin*exp(PExp_titin*abs(SL - SL_rest) - 1)*sign(SL - SL_rest)^2);

PF_collagen = heav(SL-SL_collagen)*PCon_collagen*exp(PExp_collagen*(SL-SL_collagen)-1);
dPF_collagen = heav(SL-SL_collagen)*PCon_collagen*PExp_collagen*exp(PExp_collagen*(SL-SL_collagen)-1);
   
PF = PF_titin + 1*PF_collagen; %Trabeculae
dPF = dPF_titin + dPF_collagen;
end