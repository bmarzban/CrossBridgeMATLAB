function [ dYdT, dSL, Ftotal ] = Model_XB( t,y,TmpC,MgATP,Pi,MgADP,dSL)
% Written by: Shivendra Tewari
% E-mail: TewariSG@gmail.com
% This code simulates sarcomere quick-release experiments based on the
% reduced 4-state XB model presented in "Dynamics of cross-bridge cycling, 
% ATP hydrolysis, force generation, and deformation in cardiac muscle". For
% model equations please refer to the manuscript.
% Disclaimer: This code is free to use, edit, reproduce as long as the
% source is cited.
%            P <--> 1 (t,s)
%            |      |
%      (t,s) 3 <--> 2 (t,s)

%% Constants and parameters
% Estimated parameters from Pi and ATP data (Average of N=21 GA runs)
par = [4.5397e+02   1.2521e+02   4.1169e+01   1.7553e+01   1.5928e+02   1.5372e+00   8.7750e+01   1.5137e+01   1.0060e+01   5.0247e+01   9.9383e-03   4.0067e+00   7.2899e+02   5.0129e-01    1.1370e+03   2.5464e+02   1.9066e+04   5.9698e-01];
% Estimated Q10s from individual mouse
Q10s = [1.4382e+00   3.0265e+00   1.0717e+00   1.3403e+00   1.4782e+00   1.4413e+00];

alpha1 = 1*par(8);
alpha2 = 1*par(9);
alpha3 = 1*par(10);
s3 = par(11);
K_Pi = par(12);
K_T = par(18); 

K_T = 1; 

K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T); 
g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D); 
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi); f2 = 1/(1 + Pi/K_Pi); 

kf = par(1)*Q10s(1)^((TmpC-17)/10);
kb = par(2)*f1*Q10s(1)^((TmpC-17)/10);
k1 = par(3)*f2*Q10s(1)^((TmpC-17)/10);%
k_1 = par(4)*Q10s(1)^((TmpC-17)/10);%
k2 = par(5)*1*Q10s(2)^((TmpC-17)/10);
k_2 = par(6)*1*Q10s(1)^((TmpC-17)/10)*g1;
k3 = par(7)*Q10s(1)^((TmpC-17)/10)*g2;%

kf = kf.*(1/(1+(abs(dSL)/0.4).^2));

% k3 = k3 ;

kpe2 = par(13)*Q10s(3)^((TmpC-17)/10);
eta = par(14)*Q10s(3)^((TmpC-17)/10);
kstiff1 = par(15)*Q10s(4)^((TmpC-17)/10); 
kpe1 = par(16)*Q10s(5)^((TmpC-17)/10);
kstiff2 = par(17)*Q10s(6)^((TmpC-17)/10); 

SL_max = 2.4; 
SL_min = 1.4;
SL_rest = 1.9;  % (um)

%% State Variables
P1o = y(1);
P1i = y(2);
P1w = y(3);

P2o = y(4);
P2i = y(5);
P2w = y(6);

P3o = y(7);
P3i = y(8);
P3w = y(9);
Pu = 1 - P1o - P2o - P3o;

SL = y(10);
%% Stretch-sensitive rates
f_alpha1o = (P1o - alpha1*P1i + 0.5*(alpha1*alpha1)*P1w);
f_alpha1i = (P1i - alpha1*P1w);

alpha0 = 1*alpha1;
f_alpha0o = (P2o + alpha0*P2i + 0.5*alpha0*alpha0*P2w);
f_alpha0i = (P2i + alpha0*P2w);

f_alpha2o = (P2o - alpha2*P2i + 0.5*(alpha2*alpha2)*P2w);
f_alpha2i = (P2i - alpha2*P2w);

alpha2b = 0; 
f_alphao = (P3o + alpha2b*P3i + 0.5*(alpha2b*alpha2b)*P3w);
f_alphai = (P3i + alpha2b*P3w);

f_alpha3o = (P3o + alpha3*(s3*s3*P3o + 2*s3*P3i + P3w));
f_alpha3i = (P3i + alpha3*(s3*s3*P3i + 2*s3*P3w));

%% Compute Active & Passive Force
% Active Force
dr = 0.01; % Power-stroke Size; Units: um
B_process = kstiff2*dr*P3o;   % Force due to XB cycling
C_process = kstiff1*(P2i+P3i);% Force due to stretching of XBs
F_active = (B_process + C_process);
    
% Non-linear Passive force; Adopted from Rice etal (Biophys J. 2008 Sep;95(5):2368-90)
[F_passive,dFpassive] = passiveForces(SL,SL_rest,kpe1);


% f_myofibril = 0.45; % Percent Myofibril in Muscle from Palmer etal (Mol Cell Biochem. 2004 Aug;263(1-2):73-80)

f_myofibril = 1; % Percent Myofibril in Muscle from Palmer etal (Mol Cell Biochem. 2004 Aug;263(1-2):73-80)
Ftotal = f_myofibril*(F_active + F_passive);%

%% XB ODEs
% dSL = ((intf/eta) - dfxb/kpe2)*heav(SL-SL_min)*heav(SL_max-SL)/den;
% dsL1=heav(t - 0.15)*dSL *heav(SL-SL_min)*heav(SL_max-SL);
dSL1 =0;
dP1o = kf*Pu   - kb*P1o - k1*f_alpha1o + k_1*f_alpha0o;
dP1i = 1*dSL*P1o - kb*P1i - k1*f_alpha1i + k_1*f_alpha0i;
dP1w = 2*dSL*P1i - kb*P1w - k1*P1w + k_1*P2w;

dP2o =         - k_1*f_alpha0o - k2*f_alpha2o + k_2*f_alphao + k1*f_alpha1o;
dP2i = 1*dSL*P2o - k_1*f_alpha0i - k2*f_alpha2i + k_2*f_alphai + k1*f_alpha1i;
dP2w = 2*dSL*P2i - k_1*P2w       - k2*P2w + k_2*P3w + k1*P1w;

dP3o =         + k2*f_alpha2o - k_2*f_alphao - k3*f_alpha3o;
dP3i = 1*dSL*P3o + k2*f_alpha2i - k_2*f_alphai - k3*f_alpha3i;
dP3w = 2*dSL*P3i + k2*P2w       - k_2*P3w - k3*P3w;
    
dYdT = [dP1o; dP1i; dP1w; dP2o; dP2i; dP2w; dP3o; dP3i; dP3w; dSL1];
end