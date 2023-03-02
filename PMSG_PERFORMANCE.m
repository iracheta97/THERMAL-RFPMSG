function[GEN]= PMSG_PERFORMANCE(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
wm=GEN.wr.Data;
Emax=GEN.Emax.Data;             % Volts
%--------------------------------------------------------------------------
% Topological Constraints:
%--------------------------------------------------------------------------
Nm=GEN.CONSTRAINTS.Nm;
Nph=GEN.CONSTRAINTS.Nph;
Ns=GEN.CONSTRAINTS.Ns;
%--------------------------------------------------------------------------
% Geometry:
%--------------------------------------------------------------------------
Rso=GEN.GEOMETRY.Rso;
L=GEN.GEOMETRY.L;
d3=PMSG.GEOMETRY.d3;
wsb=PMSG.GEOMETRY.wsb;
%--------------------------------------------------------------------------
% Steel:
%--------------------------------------------------------------------------
kh=GEN.STEEL.Kh;                % Hysteresis Loss Coefficient
ke=GEN.STEEL.Ke;                % Excess core loss coefficient
kc=GEN.STEEL.Ke;                % Eddy-current core loss coefficient
Vst=GEN.STEEL.Vst;              % Steel Volume
Bmax=GEN.STEEL.Bmax;            % Maximum Flux Density
% -------------------------------------------------------------------------
% ELECTRICAL PARAMETERS:
%--------------------------------------------------------------------------
Emax= GEN.ELECTRICAL.Emax;
Iph= GEN.ELECTRICAL.Iph;
Rph= GEN.ELECTRICAL.Rph;
%--------------------------------------------------------------------------
% Pout:
Pout=Nph*Emax*Iph;          % Watts
GEN.PERFORMANCE.Pout=Pout;

% Ohmic Losses:
Pr= Nph*Iph^2*Rph;
GEN.PERFORMANCE.Pr=Pr;      % Watts

% Core-Loss
fe= GEN.fe.Data;                                                % Hz
we= 2*pi*fe;
Pcl= (kh*fe*Bmax^2+ke*fe^2*Bmax^2+ke*fe^1.5*Bmax^1.5)*Vst;      % Watts
GEN.PERFORMANCE.Pcl=Pcl;

% Watts

% Stray-Losses:
Ps=0.01*Pout;
GEN.PERFORMANCE.Ps=Ps;                                          % Watts

% Input Power:
Pin=Pout+Ps+Pcl+Pr;                                             % Watts
GEN.PERFORMANCE.Pin=Pin;                                        % Watts

% Efficiency:
eta=100*Pout/Pin;
GEN.PERFORMANCE.eta=eta;                                        % Watts

% Heat Density:
qr=Pr/(L*(2*d3+wsb)*Ns);                                        % W/m2
GEN.PERFORMANCE.qr=qr;                                          % W/m2

% Maximum Heat Density:
qrt=(Pr+Pcl)/(2*pi*Rso*L);                                        % W/m2
GEN.PERFORMANCE.qrt=qrt;                                          % W/m2
