function[GEN]= PMSG_FLUX(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
%--------------------------------------------------------------------------
% Topological Constraints:
%--------------------------------------------------------------------------
Nm=GEN.CONSTRAINTS.Nm;
Nsm=GEN.CONSTRAINTS.Nsm;
Ns=GEN.CONSTRAINTS.Ns;
%--------------------------------------------------------------------------
% Geometry:
%--------------------------------------------------------------------------
g= GEN.GEOMETRY.g;
Rso= GEN.GEOMETRY.Rso;
Rsi= GEN.GEOMETRY.Rsi;
Rro= GEN.GEOMETRY.Rro;
lm= GEN.GEOMETRY.lm;
Tp=GEN.GEOMETRY.Tp;
Ts=GEN.GEOMETRY.Ts;
L=GEN.GEOMETRY.L;
am=GEN.GEOMETRY.am;
Cphi=GEN.GEOMETRY.Cphi;
Tetas=GEN.GEOMETRY.Tetas;
%--------------------------------------------------------------------------
% Magnet:
%--------------------------------------------------------------------------
miur= GEN.MAGNET.miuR;
kml= GEN.MAGNET.kml;
kc= GEN.MAGNET.kc;
PC= GEN.MAGNET.PC;
Br= GEN.MAGNET.Br;
%--------------------------------------------------------------------------
% Steel:
%--------------------------------------------------------------------------
Bmax= GEN.STEEL.Bmax;
kst= GEN.STEEL.Kst;
%--------------------------------------------------------------------------
% Shoe:
%--------------------------------------------------------------------------
as= GEN.SHOE.as;
asd= GEN.SHOE.asd;
ws= GEN.SHOE.ws;
%--------------------------------------------------------------------------
% Thermal:
%--------------------------------------------------------------------------
Txm=PMSG.THERMAL.Tm;
a=PMSG.MAGNET.a;
% -------------------------------------------------------------------------
% FLUX CALCULATIONS:
%--------------------------------------------------------------------------
% Temperature correction (20 to 120 °C):
% Br=Br*(1-a*(Txm-20)/100);                   % T
% GEN.MAGNET.Br=Br;
% Bg:
CTE= 1+ miur*kml*kc/PC;
Bg=Cphi*inv(CTE)*Br;
GEN.FLUX.Bg=Bg;            % T
% Bbi:
GEN.FLUX.Bbi=0.5*Bg;       % T
% Phig:
Ag=Tp*L*(1+am)/2;          % m2
GEN.GEOMETRY.Ag=Ag;        % m2
% Air Gap Volume:
Vg=Ag*g;                   % m3
GEN.GEOMETRY.Vg=Vg;        % m3
% Air-Gap Flux
Phig=Bg*Ag;                % Wb
GEN.FLUX.Phig=Phig;        % Wb
% Phibi:
GEN.FLUX.Phibi=0.5*Phig;   % Wb
% Stator Dimensions:
wbi= Phig/(2*Bmax*kst*L);  % 
GEN.GEOMETRY.wbi=wbi;      % m
% Slot-Bottom Radious:
Rsb= Rso-wbi;              % m
GEN.GEOMETRY.Rsb=Rsb;      % m
% Inside Rotor Radius:
Rri= Rro-lm-wbi;           % m
GEN.GEOMETRY.Rri=Rri;      % m
% Slot Depth:
ds= Rsb-Rsi;               % m
GEN.GEOMETRY.ds=ds;        % m

% Toot width:
wtb=2*wbi/Nsm;             % m
GEN.GEOMETRY.wtb=wtb;      % m

% Tooth width at Rsi
wt=Ts-ws;
GEN.GEOMETRY.wt=wt;       % m

% Slot width at Rsb:
wsb= Rsb*Tetas-wtb;
GEN.GEOMETRY.wsb=wsb;     % m

% Slot Fraction:
fwtb=wtb/(wtb+wsb);             %
fwsb= 1-fwtb;
wsi=(Rsi+asd*wtb)*Tetas-wtb;    % m
%wsi=(Rsi+asd*wtb)*Tetas*fwsb;  % m
GEN.GEOMETRY.wsi=wsi;           % m
%asd= wsi/(wsi+wtb)

% d3:
d3= ds-asd*wtb;
GEN.GEOMETRY.d3=d3;             % m

% d12
d12=asd*wtb;                    % m
d1=0.5*d12;
d2=d1;
GEN.GEOMETRY.d1=d1;             % m
GEN.GEOMETRY.d2=d2;             % m

% As: (Slot are by considering only d3)
As=d3*(Tetas*(Rsb-d3/2)-wtb);   %m2
GEN.GEOMETRY.As=As;             % m

% Slot pitch (Ts12): at Rsi+d1+d2:
Ts12=2*pi*(Rsi+d1+d2)/Ns;
GEN.GEOMETRY.Ts12=Ts12;         % m

% Tooth width at Rsi+d1+d1:
wti=Ts12-wsi;
GEN.GEOMETRY.wti=wti;           % m

% Temperature correction (20 to 120 °C):
% Brx
Brx=Br*(1+a*(Txm-20)/100);                   % T
GEN.MAGNET.Brx=Brx;
% Bg:
Bgx=Cphi*inv(CTE)*Brx;                       % T (Dependes from temperatura)
Bg=Cphi*inv(CTE)*Br;                         % T
GEN.FLUX.Bgx=Bgx;
GEN.FLUX.Bg=Bg;
% Shoe stator flux (Bts):
Bts=Bgx/(1-as);               % T
GEN.FLUX.Bts=Bts;             % T
% Tooth stator flux (Bts):
Btb=Bgx/(1-wtb/Ts);           % T
GEN.FLUX.Btb=Btb;             % T