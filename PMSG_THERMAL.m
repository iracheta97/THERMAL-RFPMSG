function[GEN]= PMSG_THERMAL(PMSG)
% Help: This function solves the thermal model of a RF/PMSG
% Created by: Dr. Reynaldo Iracheta Cortez
% Date: August 3, 2020
% Inputs:
% PMSG: Structure Data of the PMSG
% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
% Frequency:
f=PMSG.fe.Data;
% Geometry:
g=GEN.GEOMETRY.g;               % m
L=GEN.GEOMETRY.L;               % m
Rso=GEN.GEOMETRY.Rso;           % m
Rsi=GEN.GEOMETRY.Rsi;           % m
Rsb=GEN.GEOMETRY.Rsb;           % m
Rro=GEN.GEOMETRY.Rro;           % m
Rri=GEN.GEOMETRY.Rri;           % m
wtb=GEN.GEOMETRY.wtb;           % m
wt=GEN.GEOMETRY.wt;             % m
ws=GEN.SHOE.ws;                 % m
wbi=GEN.GEOMETRY.wbi;           % m
wsb=GEN.GEOMETRY.wsb;           % m
ds=GEN.GEOMETRY.ds;             % m
d3=GEN.GEOMETRY.d3;             % m
d2=GEN.GEOMETRY.d1;             % m
d1=GEN.GEOMETRY.d1;             % m
Ts=GEN.GEOMETRY.Ts;             % m
Tp=GEN.GEOMETRY.Tp;             % m
am=GEN.GEOMETRY.am;             % m
lm=GEN.GEOMETRY.lm;             % m
As=GEN.GEOMETRY.As;             % m
% Coil-Pitch Fraction:
acp=GEN.GEOMETRY.acp;              %
% Topological Constraints:
Nm =GEN.CONSTRAINTS.Nm;
Ns =GEN.CONSTRAINTS.Ns;
Nsm=GEN.CONSTRAINTS.Nsm;
Nsp=GEN.CONSTRAINTS.Nsp;
% Factors:
Kst=GEN.STEEL.Kst;
% Thermal Radiation constants:
a1= GEN.THERMAL.asy;                % W/(°K-m2)
a2= GEN.THERMAL.ag;                 % W/(°K-m2)
a3= GEN.THERMAL.aes;                % W/(°K-m2)
a4= GEN.THERMAL.aew;                % W/(°K-m2)
a5= GEN.THERMAL.aew;                % W/(°K-m2)
aframe= GEN.THERMAL.afr;            % W/(°K-m2)
% Thermal Conduction constants:
lamair=GEN.THERMAL.lamair;          % W/(°K-m)
lamfe=GEN.THERMAL.lamfe;            % W/(°K-m)
lami=GEN.THERMAL.lami;              % W/(°K-m)
lamcoil=GEN.THERMAL.lamcoil;        % W/(°K-m)
lamcu=GEN.THERMAL.lamcu;            % W/(°K-m)
lamgpr=GEN.THERMAL.lamgrp;          % W/(°K-m)
lamm=GEN.THERMAL.lamm;              % W/(°K-m)
lamglue=GEN.THERMAL.lamglue;        % W/(°K-m)
lamf=GEN.THERMAL.lamf;              % W/(°K-m)
% Cooling air constants:
rhoc=GEN.THERMAL.rhoc;              % kg/m3
kthc=GEN.THERMAL.kthc;              % J/(kg-°K)
vair=GEN.THERMAL.vair;              % m/s
hair=GEN.THERMAL.hair;              % m
% Frame:
hframe=GEN.THERMAL.hframe;          % m
% Flux densities:
Bysmax=GEN.STEEL.Bmax;              % T
Byrmax=Bysmax;                      % T
Bg=GEN.FLUX.Bg;                     % T
Bbi=GEN.FLUX.Bbi;                   % T
Bts=GEN.FLUX.Bts;                   % T
Btb=GEN.FLUX.Btb;                   % T
% Air temperature:
Tair=PMSG.THERMAL.Tair;             % °C
% Contact surfaces:
Tso=2*pi*(Rso)/Ns;
wso=Tso-wtb;
bd=wtb;
bs=wso;
% Winding:
kcp=GEN.WINDING.kcp;               % dimensionless
le2=GEN.WINDING.le2;               % m
le1=GEN.WINDING.le1;               % m
% Thickness (For thermal analysis):
lgrp=GEN.THERMAL.lgrp;               % m
lmglue=GEN.THERMAL.lmglue;           % m
lci=GEN.THERMAL.lci;                 % m
% Thermal resistances per pole:
ti=lci;                                    % m
bcu=wsb-2*ti;                              % m
GEN.WINDING.bcu=bcu;                       % m
te=wtb+2*ti;
GEN.WINDING.te=te;                         % m
dcu=d3-2*ti;                               % m
% Length of the diagonal portion of the end winding:
le3=(acp*Tp/2)*Ts/(sqrt(Ts^2-(bcu+te)^2));   % m
% Length of the end-winding side:
lend1=2*(le2+abs(le3));                     % m
lb=lend1;
% Electrical parameters:
Nph=GEN.CONSTRAINTS.Nph;                    % ¨hases
Rph=GEN.ELECTRICAL.Rph;                     % Ohms
Iphx=GEN.ELECTRICAL.Iphx;                   % A
Is=GEN.ELECTRICAL.Is;                       % A
Jc=GEN.ELECTRICAL.Jc;                       % A/m2
ns=GEN.ELECTRICAL.ns;                       % Turns
np=GEN.ELECTRICAL.np;                       % Turns
Ac=GEN.ELECTRICAL.Ac;
Vcu=GEN.WINDING.Vcu;                       % m3
rhocu75=GEN.WINDING.rhocu75;               % Ohms-m
rhocu=GEN.WINDING.rhocu;                   % Ohms-m
Tc=GEN.GEOMETRY.Tc;                        % m
% Steel:
rhobi=GEN.STEEL.rhobi;                     % kg/m3
% Magnet:
lm=PMSG.GEOMETRY.lm;                       % m
Wm=PMSG.MAGNET.Wm;                         % kg
Wst=PMSG.STEEL.Wst;                        % kg
Wro=PMSG.STEEL.Wro;                        % kg
Wcu=PMSG.WINDING.Wcu;                      % kg
Am=PMSG.GEOMETRY.Am;                       % m2
%--------------------------------------------------------------------------
% Radiation:
%--------------------------------------------------------------------------
% Frame:
Rfr=1/(aframe*pi*(hframe+2*Rso)*L);
% Stator yoke:
R1=1/(3*bd*L*a1);               % W/°K-m
R2=1/(3*bs*L*a1);               % W/°K-m
% Shoes:
%R22= (d1+d2)/(L*(Ts-ws)*a2);   % W/°K-m
R22= 1/(L*(Ts-ws)*a2);          % W/°K-m
% Rotor yoke:
R28=1/(L*Tp*a5);                % W/°K-m
% Axial: PMSG End-shield:
R29=1/(pi*Rso^2*a3);            % W/°K-m
% Axial: End windings:
R35=2/(lb*bcu*a4);              % W/°K-m
R38=2/(lb*dcu*a4);              % W/°K-m
%--------------------------------------------------------------------------
% Convection:
%--------------------------------------------------------------------------
% Temperature rise of cooling air ducts:
% Number of cooling circuits:
tair=Rso/6;                                % m
Nair=pi*Rso/tair;
% Volumetric air flow in each channel:
qvair=vair*L*hair;                         % m3/s
% Total volumetric cooling air flow:
qvc=Nair*qvair;                            % m3/s
Rfcv=1/(qvc*rhoc*kthc);                    % W/°K-m
R0=Rfr*Rfcv/(Rfr+Rfcv);                    % W/°K-m
PMSG.THERMAL.R0=R0;
%--------------------------------------------------------------------------
% Conduction:
%--------------------------------------------------------------------------
% Stator yoke:
lu=L*Kst;                                  % Effective length of laminations
Rbi=wbi/2;                                 % m
R3=0.5*(wbi+Rbi)/(lamfe*bd*lu);            % W/°K-m
R4=0.5*(wbi+0.5*Rbi)/(lamfe*bs*lu);        % W/°K-m
R5=0.5*bd/(lamfe*(wbi+Rbi)*lu);            % W/°K-m
R6=0.5*bs/(lamfe*(wbi+Rbi)*lu);            % W/°K-m
R7=-R5/3;                                  % W/°K-m
R8=-R6/3;                                  % W/°K-m
R9=-R3/3;                                  % W/°K-m
R10=-R4/3;                                 % W/°K-m
%
R7=R5/3;                                  % W/°K-m
R8=R6/3;                                  % W/°K-m
R9=R3/3;                                  % W/°K-m
R10=R4/3;                                 % W/°K-m
%%
R7=R5/3;                                  % W/°K-m
R8=R6/3;                                  % W/°K-m
R9=R3/3;                                  % W/°K-m
R10=R4/3;                                 % W/°K-m
% Conductors (Radial):
R11=ti/(L*bcu*lami);                       % W/°K-m
R12=0.5*(d3)/(lu*bd*lamfe);                % W/°K-m
R13=0.5*dcu/(lu*bcu*lamcoil);              % W/°K-m
R14=0.5*bd/(lu*d3*lamfe);                  % W/°K-m
R15=0.5*ti/(lu*dcu*lami);                  % W/°K-m
R16=0.5*bcu/(lu*dcu*lamcoil);              % W/°K-m
R17=-R14/3;                                % W/°K-m
R18=-R16/3;                                % W/°K-m
R19=-R12/3;                                % W/°K-m
R20=-R13/3;                                % W/°K-m
%%
R17=R14/3;                                % W/°K-m
R18=R16/3;                                % W/°K-m
R19=R12/3;                                % W/°K-m
R20=R13/3;                                % W/°K-m
% Shoes:
R21= 0.5*(d1+d2)/(lu*(0.5*wt+0.5*bd)*lamfe);   % W/°K-m
% Magnet:
bm=Tp*am;                                   % m
R23=1/(L*bm*a2)+ lgrp/(L*bm*lamgpr);        % W/°K-m
R24=0.5*lm/(L*bm*lamm);                     % W/°K-m
%R25=-R24/3;                                % W/°K-m
R25=R24/3;                                  % W/°K-m
R26=lmglue/(L*bm*lamglue);                  % W/°K-m
% Rotor yoke:
hyr=wbi;                                    % m
R27=hyr/(L*Tp*lamfe);                       % W/°K-m
% Conductors (Axial):
R31=0.5*L/(dcu*bcu*kcp*lamcu);             % W/°K-m
R30=-R31/3;                                % W/°K-m
R30=R31/3;                                 % W/°K-m
R32=0.5*lb/(dcu*bcu*kcp*lamcu);            % W/°K-m
R33=-R32/3;
R33=R32/3;
R36=0.5*dcu/(lb*bcu*lamcoil);              % W/°K-m
R34=-R36/3;
R34=R36/3;
R39=0.5*bcu/(lb*dcu*lamcoil);              % W/°K-m
R37=-R39/3;                                % W/°K-m
R37=R39/3;                                 % W/°K-m
% Air-Gap Thermal Resistance:
R40=0.5*g/(lamair*pi*L*(Rsi+Rro)/2);       % W/°K-m
% Full model:
RF1=[R0;R1;R2;R3;R4;R5;R6;R7;R8;R9;R10;R11;R12;R13;R14;R15;R16;R17;R18;R19;R20;...
    R21;R22;R23;R24;R25;R26;R27;R28;R29;R30;R31;R32;R33;R34;R35;R36;R37;R38;R39;R40];
PMSG.THERMAL.RF1=RF1;                       % W/°K-m
%--------------------------------------------------------------------------
% Full simplified model:
%--------------------------------------------------------------------------
R50=R0;
R51=(R1+R3)/Ns;
R52=(R2+R4)/Ns;
R53=R9/Ns;
R54=(R7+R8+0.5*(R5+R6))/Ns;
R55=R10/Ns;
R56=(R3+R12)/Ns;
R57=(R4+R11+R13)/Ns;
R58=R19/Ns;
R59=(R17+R18)/Ns+0.5*(R14+R15+R16)/Ns;
R60=R20/Ns;
R61=(R12+R21)/Ns;
R62=(R13+R11)/Ns;
R63= (R21+R22)/Ns;                          % W/°K-m
R64=R40;                                    % W/°K-m
R65=R40;                                    % W/°K-m
R66=(R24+R23)/Nm;                           % W/°K-m
R67=(R25)/Nm;                               % W/°K-m
R68=(R24+R26)/Nm;                           % W/°K-m
R69= R27/Nm;                                % W/°K-m
R70=R30/Ns;
R71=(R33+R32+R32)/Ns;                       % W/°K-m
R72=R68;                                    % W/°K-m
R73a= (0.5*(R35+R36)+R34)/Ns;               % W/°K-m
R73b= (0.5*(R35+R36)+R34)/Ns;               % W/°K-m
R73= R73a*R73b/(R73a+R73b);                 % W/°K-m
R74=R73;                                    % W/°K-m
R75=2*R29;
R76=R75;                                    % W/°K-m
RF2=[R50;R51;R52;R53;R54;R55;R56;R57;R58;R59;R60;R61;R62;R63;R64;R65;R66;R67;R68;R69;R70;...
    R71;R72;R73;R74;R75];
PMSG.THERMAL.RF2=RF2;                       % W/°K-m
%--------------------------------------------------------------------------
% Matrix of thermal resistances:
%--------------------------------------------------------------------------
RRF=[R50 1 0;
    R51 2 1;
    R52 5 1;
    R53 3 2;
    R54 4 3;
    R55 5 4;
    R56 6 2;
    R56 11 10;
    R57 9 5;
    R58 7 6;
    R59 9 8;
    R60 9 8;
    R61 10 6;
    R62 11 9;
    R63 12 10;
    R64 12 11;
    R65 13 11;
    R66 14 13;
    R67 20 14;
    R68 15 14;
    R69 16 15;
    R70 19 8;
    R71 19 17;
    R72 19 18;
    R73 17 16;
    R74 18 16;
    R75 16 0;
    R76 16 0];
% Resistances
R=RRF(:,1);
%POS=1+[0 -1;1 0;2 0;1 2; 1 3;2 5;9 10; 4 5; 9 10;3 4;4 5;3 9;5 10;9 11;11 8;8 -1;4 6;10 6;7 8;6 7];
POS=RRF(:,2:3);
GEN.THERMAL.RRF=RRF;                         % Full model
nodes=max(max(POS));
GEN.THERMAL.NODEF=nodes;                     % Full model 
%--------------------------------------------------------------------------
% Power Losses:
%--------------------------------------------------------------------------
% Copper losses:
Pcu=Nph*Rph*(Iphx)^2;                       % W
PMSG.THERMAL.Pcu=Pcu;                       % W
% Core losses:
% Stator:
Vys=(pi*(Rso^2-Rsb^2)-Ns*pi*Rbi^2)*L*Kst;           % m3
Wys=rhobi*Vys;                                      % kg
GEN.THERMAL.Wys=Wys;                                % kg
% Tooth:
Vt=(pi*((Rsb)^2-(Rsi+d1+d2)^2)-Ns*As)*L*Kst;        % m3
GEN.THERMAL.Vt=Vt;                                  % m3
Wt=rhobi*Vt;                                        % kg
GEN.THERMAL.Wt=Wt;                                  % kg
% Vt2=Ns*wt*d3*L*Kst;                               % m3
% GEN.THERMAL.Vt2=Vt2;                              % m3
% Wt2=rhobi*Vt2;                                    % kg
% GEN.THERMAL.Wt2=Wt2;
% Shoe:
Vsh=Ns*(d1+d2)*(0.5*wt+0.5*bd)*L*Kst;               % m3
Wsh=rhobi*Vsh;
GEN.THERMAL.Wsh=Wsh;
% Hysteresis losses:
% Specific hysteresis and eddy current losses (50 Hz and 1.5 T):
phy=2.04;                                  % W/kg
pft=0.76;                                  % W/kg
% Stator: Empirical factors (Hysteresis + Eddy currents):
khy=2;
kft=1.8;
Bysmax=GEN.STEEL.Bmax;              % T
Byrmax=Bysmax;                      % T
Bgx=GEN.FLUX.Bgx;                   % T
Bbi=GEN.FLUX.Bbi;                   % T
%Bbi=Bysmax;
Bts=GEN.FLUX.Bts;                   % T
%Bts=Bysmax;
Btb=GEN.FLUX.Btb;                   % T
%Btb=Bysmax;
phys=khy*Wys*phy*(f/50)*(Bbi/1.5)^2;       % W
GEN.THERMAL.phys=phys;
pfts=kft*Wys*pft*(f/50)^2*(Bbi/1.5)^2;     % W
GEN.THERMAL.pfts=pfts;
% Teeth: Empirical factors (Hysteresis + Eddy currents):
khyd=1.2;
kftd=2.5;
phyd=khyd*Wt*phy*(f/50)*(Btb/1.5)^2;       % W
GEN.THERMAL.phyd=phyd;
pftd=kftd*Wt*pft*(f/50)^2*(Btb/1.5)^2;     % W
GEN.THERMAL.pftd=pftd;
% Shoe Losses:
phyz=khyd*Wsh*phy*(f/50)*(Bts/1.5)^2;       % W
GEN.THERMAL.phys=phyz;
pftz=kftd*Wsh*pft*(f/50)^2*(Bts/1.5)^2;     % W
GEN.THERMAL.pfts=pfts;
% Magnet:
pfm=300;                                  % W/m2
pftm=pfm*Nm*Am*lm;                        % W
GEN.THERMAL.pftm=pftm;                    % W
% Aditional losses:
% They are considered 20% of the core losses:
pad=0.2*(phys+pfts+phyd+pftd);             % W
GEN.THERMAL.pad=pad;                       % W
% Friction and winding losses:
% 0.5% of PN
PN=GEN.Pout.Data;                          % W
pfw=0.005*PN;                              % W
GEN.THERMAL.pfw=pfw;                       % W
% Total losses:
GEN.THERMAL.Plossm=pad/0.2+pfw+pad+pftm;   % W
%--------------------------------------------------------------------------
% Power Losses (Simplified Model):
%--------------------------------------------------------------------------
% Thermal resistances:
% Heat sources (Full Model):
%Tair=30;                                  %
% Stator yoke:
Pysa=(wtb/(wtb+wsb))*(phys+pfts);          % W
Pysb=(wsb/(wtb+wsb))*(phys+pfts);          % W
% Teeth:
Pth=(phyd+pftd);                           % W
% Shoes:
Psh=(phyz+pftz);                           % W
% Coper losses in L:
Pcua=Pcu*(L/(L+lb));                       % W
% End Winding:
Pew=Pcu*(lb/(L+lb));                        % W
% Magnet:
Pm=pftm;                                    % W
PS1=0+Tair/R50;                             % W
PS2=Pysa;                                   % W
PS3=Pysb;                                   % W
PS4=(phyd+pftd);                            % W
PS5=Pcu*(L/(L+lb));                         % W
PS6=0;                                      % W
PS7=(phyz+pftz);                            % W 
PS8=0;                                      % W
PS9=pftm;                                   % W 
PS10=0;                                     % W
PS11=Pcu*(lb/(L+lb));                       % W
PS12=0;                                     % W
%--------------------------------------------------------------------------
% Power Losses (Full model):
%--------------------------------------------------------------------------
%P0=0;                                      % W
P1=0+Tair/R50;                             % W
P2=0;                                      % W
P3=Pysa;                                   % W
P4=Pysb;                                   % W
P5=0;                                      % W
P6=0;                                      % W
P7=Pth;                                    % W
P8=Pcua;                                   % W
P9=0;                                      % W
P10=Psh;                                   % W
P11=0;                                     % W
P12=0;                                     % W
P13=0;                                     % W
P14=0;                                     % W
P15=0;                                     % W
P16=0+Tair/R29;                            % W
P17=Pew;                                   % W
P18=Pew;                                   % W
P19=0;                                     % W
P20=Pm;                                    % W
PLOSS=[P1;P2;P3;P4;P5;P6;P7;P8;P9;P10;P11;P12;P13;P14;P15;P16;P17;P18;P19;P20]; % W
GEN.THERMAL.PF=PLOSS;                      % W
%--------------------------------------------------------------------------
% Thermal conductance matrix: (FULL MODEL)
%--------------------------------------------------------------------------
G= zeros(nodes,nodes);
% 
RR=isempty(R);
if RR==0
% Resistencias del sistema:
NR= length(R);
% Calculo de las conductancias del sistema:
if NR>0
GR=1./R;
% Verificación de datos:
NX= length(POS(:,1));
NY= length(POS(:,2));
if NR==NX & NR==NY
    % Inicialización de las historias del sistema:
    for n=1:NR
        X= POS(n,1);
        Y= POS(n,2);
            if  X==0 | Y==0
                if X>Y
                G(X,X)= GR(n)+ G(X,X);
                else
                G(Y,Y)= GR(n)+ G(Y,Y);  
                end
            else
                G([X Y],[X Y])=[1 -1; -1 1]*GR(n)+G([X Y],[X Y]);
            end         
    end
                 
else
    error('ERROR EN LAS CONEXIONES DE LAS RESISTENCIAS');
end

    else
    % No hay resistencias en el sistema:
    R=[];
end
    else
    % No hay resistencias en el sistema:
    R=[];
end
% Matrix of thermal conductances:
GEN.THERMAL.GF=G;
% Temperatures:
TF=inv(G)*PLOSS;
GEN.THERMAL.TF=TF;                    % °C
%--------------------------------------------------------------------------
% Matrix of thermal conductances (Simplified model):
%--------------------------------------------------------------------------
% Thermal resistances:
R50=R0;
R51=(R1+R3)/Ns;
R52=(R2+R4)/Ns;
R53=(R7+R8+R9+R10+0.5*(R5+R6))/Ns;
R54=(R3+R12)/Ns;
R55=(R4+R11+R13)/Ns;
R56=(R19+R17+R18)/Ns+0.5*(R14+R15+R16)/Ns;
R57=R20/Ns;
R58=(R12+R21)/Ns;
R59=(R13+R11)/Ns;
% Modificación:
R60=(R21+R22)/Ns+R40;
%R60=(R12+R21+R22+R40)/Ns+(R23+R24+R40)/(Nm);
R61=(R23+R25)/Nm+R40;
R62=(R24+R26+R27+R28)/Nm;
R63=(R30+0.5*(R31+R32))/Ns;
R64=0.5*R33/Ns;
R65a=0.5*(R34+0.5*(R36+R35))/Ns;
R65b=0.5*(R37+0.5*(R39+R38))/Ns;
R65=R65a*R65b/(R65a+R65b);
R66=R29;
RR=[R50 1 0;
    R51 2 1;
    R52 3 1;
    R53 3 2;
    R54 4 2;
    R55 6 3;
    %R56 5 4;
    R56 5 4;
    R57 6 5;
    R58 7 4;
    R59 6 8;
    R60 7 8;
    R61 9 8;
    R62 9 10;
    R63 12 5;
    R64 11 12;
    R65 10 11;
    R66 10 0];
%R=[R50 R51 R52 R53 R54 R55 R56 R56 R57 R58 R59 R60 R61 R62 R63 R63 R64 R65];
R=RR(:,1);
%POS=1+[0 -1;1 0;2 0;1 2; 1 3;2 5;9 10; 4 5; 9 10;3 4;4 5;3 9;5 10;9 11;11 8;8 -1;4 6;10 6;7 8;6 7];
POS=RR(:,2:3);
GEN.THERMAL.RR=RR;                         % 
nodes=max(max(POS));
GEN.THERMAL.NODES=nodes;
%--------------------------------------------------------------------------
% Heat Sources:
%--------------------------------------------------------------------------
%Tair=30;                                   %

P1=0+Tair/R50;                             % W
P2=Pysa;                                   % W
P3=Pysb;                                   % W
P4=(phyd+pftd);                            % W
P5=Pcu*(L/(L+lb));                         % W
P6=0;                                      % W
P7=(phyz+pftz);                            % W 
P8=0;                                      % W
P9=pftm;                                   % W 
P10=0+Tair/R66;                            % W
P11=2*Pcu*(lb/(L+lb));                       % W
P12=0;                                     % W

P=[P1;P2;P3;P4;P5;P6;P7;P8;P9;P10;P11;P12]; % W
% Vector of losses:
GEN.THERMAL.PS=P;                           % W
%--------------------------------------------------------------------------
% Thermal conductance matrix:
%--------------------------------------------------------------------------
G= zeros(nodes,nodes);
% 
RR=isempty(R);
if RR==0
% Resistencias del sistema:
NR= length(R);
% Calculo de las conductancias del sistema:
if NR>0
GR=1./R;
% Verificación de datos:
NX= length(POS(:,1));
NY= length(POS(:,2));
if NR==NX & NR==NY
    % Inicialización de las historias del sistema:
    for n=1:NR
        X= POS(n,1);
        Y= POS(n,2);
            if  X==0 | Y==0
                if X>Y
                G(X,X)= GR(n)+ G(X,X);
                else
                G(Y,Y)= GR(n)+ G(Y,Y);  
                end
            else
                G([X Y],[X Y])=[1 -1; -1 1]*GR(n)+G([X Y],[X Y]);
            end         
    end
                 
else
    error('ERROR EN LAS CONEXIONES DE LAS RESISTENCIAS')
end

    else
    % No hay resistencias en el sistema:
    R=[];
end
    else
    % No hay resistencias en el sistema:
    R=[];
end
% Matrix of thermal conductances:
GEN.THERMAL.GS=G;
% Temperatures:
TS=inv(G)*P;
GEN.THERMAL.TS=TS;                    % °C
%--------------------------------------------------------------------------
% Specific temperatures:
%--------------------------------------------------------------------------
% Air temperature:
Tair=30;                                       % °C
% Copper Temperature:                          % °C
Tcu=TF(19,:);                                  % °C
GEN.THERMAL.Tcu=Tcu;                            % °C
% Magnet Temperature:
Tm=TF(14,:);                                   % °C
GEN.THERMAL.Tm=Tm;                             % °C
% Gap Temperature:
Tgap=TF(13,:);                                 % °C
GEN.THERMAL.Tgap=Tgap;                         % °C
% Stator Temperature:
Tst=TF(6,:);                                   % °C
GEN.THERMAL.Tst=Tst;                           % °C
