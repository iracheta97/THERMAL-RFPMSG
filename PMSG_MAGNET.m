function[GEN]= PMSG_MAGNET(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
%--------------------------------------------------------------------------
% Topological Constraints:
%--------------------------------------------------------------------------
Nm=GEN.CONSTRAINTS.Nm;
%--------------------------------------------------------------------------
% Geometry:
%--------------------------------------------------------------------------
g= GEN.GEOMETRY.g;
Tp=GEN.GEOMETRY.Tp;
Ts=GEN.GEOMETRY.Ts;
L=GEN.GEOMETRY.L;
am=GEN.GEOMETRY.am;
Cphi=GEN.GEOMETRY.Cphi;
%--------------------------------------------------------------------------
% Magnet:
%--------------------------------------------------------------------------
miur= GEN.MAGNET.miuR;
%--------------------------------------------------------------------------
% Shoe:
%--------------------------------------------------------------------------
as= GEN.SHOE.as;
% -------------------------------------------------------------------------
% MAGNET SIZING:
%--------------------------------------------------------------------------
% BHmax:
BHmax=GEN.MAGNET.BHmax*1000; % J/m3
% Volume:
Tem=GEN.Tem.Data;            % N-m
Vm= (2/3)*Tem/(BHmax*Nm);    % m3
GEN.MAGNET.Vm=Vm;            % m3
% Magnet Weigth:
rhom=GEN.MAGNET.rhom;        % kg/m3
Wm=Nm*Vm*rhom;               % kg
GEN.MAGNET.Wm=Wm;            % kg
% Magnet length:
Am=(L*am*Tp);
GEN.GEOMETRY.Am=Am;            % m2
lm= Vm/(Am);                   % m
GEN.GEOMETRY.lm=lm;            % m

% MAGNET ASPECT RATIO VERIFICATION:
MAR=lm/(am*Tp);
GEN.MAGNET.MAR=MAR;            %
MAR=lm/(am*Tp);
GEN.MAGNET.MAR=MAR;            % 

% PERMEANCE COEFFICIENT:
PC=lm/(g*Cphi);               %
GEN.MAGNET.PC=PC;             % 
% Ratio lm/g:
LMG=lm/g;                     %
GEN.MAGNET.LMG=LMG;           % 

% FACTOR KML:
kml=1+ 4*MAR/(pi*miur)*log(1+pi*g/((1-am)*Tp));
GEN.MAGNET.kml=kml;           %

% CARTER'S COEFFICIENT:
% Slot opening:
ws=Ts*as;
GEN.SHOE.ws=ws;             %
CTE= 1-inv(Ts*(5*g/ws+1)/ws);
kc=inv(CTE);
GEN.MAGNET.kc=kc;           %

