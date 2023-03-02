function[GEN]= PMSG_FACTORS(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
%--------------------------------------------------------------------------
% Topological Constraints:
%--------------------------------------------------------------------------
Nspp=GEN.CONSTRAINTS.Nspp;
%--------------------------------------------------------------------------
% Geometry:
%--------------------------------------------------------------------------
% Coil-Pitch Fraction:
acp=GEN.GEOMETRY.acp;              %
% Tetas:
Tetas=GEN.GEOMETRY.Tetas;           % Mechanical in Radians
% Tetase:
Tetase=GEN.GEOMETRY.Tetase;         % Electrical in Radians
% -------------------------------------------------------------------------
% FACTORS:
%--------------------------------------------------------------------------
% kd:
kd= sin(Nspp*Tetase/2)/(Nspp*sin(Tetase/2));
GEN.FACTORS.kd=kd;
% kp=acp:
kp=acp;
%kp=1;
% Coil pitch:
GEN.FACTORS.kp=kp;
% ks:
ks= 1-Tetase/(2*pi);
ks=1;
GEN.FACTORS.ks=ks;

