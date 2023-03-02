function[PMSG]= PMSG_PARAMETERS()
% Program for providing the INPUTS of the PMSG Design Program
% Created  by: Reynaldo Iracheta Cortez.
% Date: July 9th, 2017
% Inputs:

% Output
% PMSG: Structure Data of the PMSG

% clc;
% clear;
%--------------------------------------------------------------------------
% Plate-Data:
% Rated Ouput Power:
%--------------------------------------------------------------------------
Pout=10000;             % W
% Rated Speed:
Sr= 225;                % rpm
% Rated Speed:
wr= 2*pi*Sr/60;         % rad/s
% Rated Torque:
Tem= Pout/wr;           % N-m
% Start Torque: (Cogging Torque)
Tstart= 6.8;            % N-m
% Maximum Induced Voltage:
Emax=250;               % Volts
%Emax=250/sqrt(3);
% Frequency:
fe=60;                  % Hz
% Total Weigth:
W=[];
% Estimated Rotor Inertial Moment (kg-m2):
J=[];
% ARRAY STRUCTURE FOR RATED DATA:
PMSG.Pout=[];
PMSG.Pout = setfield(PMSG.Pout,'Data',Pout);
PMSG.Pout= setfield(PMSG.Pout,'Name','Rated Output Power (W)');
PMSG.Sr=[];
PMSG.Sr = setfield(PMSG.Sr,'Data',Sr);
PMSG.Sr = setfield(PMSG.Sr,'Name','Rated Speed (rpm)');
PMSG.wr=[];
PMSG.wr = setfield(PMSG.wr,'Data',wr);
PMSG.wr = setfield(PMSG.wr,'Name','Rated Speed (rad/s)');
PMSG.Tem=[];
PMSG.Tem = setfield(PMSG.Tem,'Data',Tem);
PMSG.Tem = setfield(PMSG.Tem,'Name','Rated Electromagnetic Torque (N-m)');
PMSG.Tstart=[];
PMSG.Tstart = setfield(PMSG.Tstart,'Data',Tstart);
PMSG.Tstart = setfield(PMSG.Tstart,'Name','Start Torque (N-m)');
PMSG.Emax=[];
PMSG.Emax = setfield(PMSG.Emax,'Data',Emax);
PMSG.Emax= setfield(PMSG.Emax,'Name','Maximum Induced Voltage (V)');
PMSG.fe=[];
PMSG.fe = setfield(PMSG.fe,'Data',fe);
PMSG.fe = setfield(PMSG.fe,'Name','Frequency (Hz)');
PMSG.W=[];
PMSG.W = setfield(PMSG.W,'Data',W);
PMSG.W = setfield(PMSG.W,'Name','Total Weigth (kg)');
PMSG.J=[];
PMSG.J = setfield(PMSG.J,'Data',J);
PMSG.J = setfield(PMSG.J,'Name','Estimated Rotor Inertial Moment (kg-m2)');
%--------------------------------------------------------------------------
% TOPOLOGICAL CONSTRAINTS:
%--------------------------------------------------------------------------
% Number of phases:
Nph=3;                    
% Number of poles:
Nm=32; 
% Number of layers:
Nlay=1; 
% Number of slots per phase:
Nsp=[];                    
% Number of slots:
Ns=90;
% Number of slots per phase per pole:
Nspp=[];
% Number of slots per poles:
Nsm=[];
% Number of turns per phase per pole:
ntpp=[];
PMSG.CONSTRAINTS=struct('Nph',Nph,'Nm',Nm,'Nlay',Nlay,'Nsp',Nsp,'Ns',Ns,'Nspp',Nspp, 'Nsm',Nsm,'ntpp',ntpp);
%--------------------------------------------------------------------------
% GEOMETRICAL PARAMETERS:
%--------------------------------------------------------------------------
% Air Gap:
g=1e-3;                    % m
% Outside Stator Radius:
Rso=0.1885;                 % m
% Inside Stator Radius:
Rsi=[];                     % m
% Outside Rotor Radius
Rro=0.1490;                 % m
% Inside Stator Radius:
Rri=[];                     % m
% Generator Axial Length
L=0.1955;                   % m
% Magnet Fraction:
am=[];                      % m
% Flux Concentration Factor:
Cphi=[];                    % m
% Pole Pitch:
Tp=[];                      % m
% Angular Pole-Pitch:
Tetap=[];                   % radians
% Slot-Pitch at Rsi:
Ts=[];                      % m
% Slot-Pitch ar Rsi+d1+d2:
Ts12=[];                    % m
% Angular Slot-Pitch at Rsi:
Tetas=[];                   % Mechanical radians
% Angular Slot-Pitch at Rsi in Electrical Radians:
Tetase=[];                  % Electrical radians
% Coil-Pitch at Rsi:
Tc=[];                     % m
% Coil-pitch Fraction:
acp=[];
% Coil-Pitch in Radians:
Tetace=[];
% Magnet Area:
Am=[];                      % m2
% Magnet length:
lm=[];                      % m
% Air-gap Area:
Ag=[];                      %m2
% Air-Gap Volume:
Vg=[];
% Wbi:
wbi=[];                     % m
% Slot-Bottom Radius:
Rsb=[];                     % m
% Slot depth:
ds=[];
% Tooth Width Bottom
wtb=[];
% Toot width at Rsi+d1+d2:
wti=[];
% Toot width at Rsi:
wt=[];
% Slot width at Rsb:
wsb=[];
% Slot width at Rsi+d1+d2:
wsi=[];
% d3:
d3=[];
% d1:
d1=[];
% d2:
d2=[];
% Slot Cross-Sectional Area:
As=[];

PMSG.GEOMETRY=struct('g',g,'Rso',Rso,'Rsi',Rsi,'Rro',Rro,'Rri',Rri,'L',L,'am',am,'Cphi',Cphi,...
                     'Tp',Tp,'Tetap',Tetap,'Ts',Ts,'Ts12',Ts12,'Tetas',Tetas,'Tetase',Tetase,'Tc',Tc,'acp',acp,...
                     'Tetace',Tetace,'Am',Am,'lm',lm,...
                      'Ag',Ag,'Vg',Vg,'wbi',wbi,'Rsb',Rsb,'ds',ds,'wtb',wtb','wt',wt','wti',wti,'wsb',wsb,...
                     'wsi',wsi,'d3',d3, 'd1', d1, 'd2', d2, 'As',As);
%--------------------------------------------------------------------------
% MAGNET:
%--------------------------------------------------------------------------
% Magnet Type:
Type='NdFeB N35';
% Relative Permittivity:
er= 1;                      
% Magnet Recoil Permeability:
miuR=1.05;                  %
% Residual Flux Density:
Br=1.21;                    % T (Wb/m2)
% Residual Flux Density at Tx (°C):
Brx=Br;                     % T(Wb/m2)
% Coercitive Force (kA/m):
Hc= 915;                    % kA/m
% Intrinsic Coercitive Force
Hci= 955;                   % kA/m
% Temperature coefficient of residual magnetization:
a=-0.12;
% Temperature coefficient of intrinsic coercitivity:
b=-0.6;
% Maximum Energy Product (kJ/m3)
BHmax= 270;                 % kJ/m3
% Maximum Operating Temperature:
Tmax= 80;                   % °C
% Curie Temperature:
Tcurie= 310;                % °C
% Magnet Density:
rhom= 7400;                 % kg/m3
% Magnet Volume:
Vm=[];                      % m3
% Magnet Weigth:
Wm=[];                      % kg
% Permeance Coefficient:
PC=[];                      % PC>=1
% Magnet Aspect Ratio
MAR=[];                     % MAR must be less than 0.25
% Ratio lm/g
LMG=[];                     % lm/g >1
% Magnet Leakage Factor:
kml=[];                     % 
% Carter's Coefficient:
kc=[];
% Carter's gap:
gc=[];
% Effective gap:
gef=[];
PMSG.MAGNET=struct('Type',Type,'miuR',miuR,'Br',Br,'Brx',Brx','Hc',Hc,'Hci',Hci,'BHmax',BHmax,...
                   'Tmax',Tmax,'Tcurie',Tcurie,'rhom',rhom,'Vm',Vm,'Wm',Wm,'PC',PC,...
                   'MAR',MAR,'LMG',LMG,'kml',kml,'kc',kc,'gc',gc,'gef',gef,'a',a,'b',b);
%--------------------------------------------------------------------------
% STEEL:
%--------------------------------------------------------------------------
% Steel Type:
Type='M36_24G';
% Thickness:
Thickness=0.35e-3;           % m
% Maximum Steel-Density:
Bmax=2.0088;
% Steel Permeability:
miur=28697;                  %
% Lamination Stacking Factor:
Kst=0.95;
% Core loss coefficient:
Kh= 211.6;
% Excess Core Loss Coefficient:
Kex= 1.51;
% Eddy-Current Core Loss Coefficient:
Ke= 1.14;
% Steel Density:
rhobi= 7850;                 % kg/m3
% Stator Volume:
Vst=[];                      % m3
% Stator Weigth:
Wst=[];                      % kg
% Rotor Volume:
Vro=[];                      % m3
% Rotor Weigth:
Wro=[];                      % kg

PMSG.STEEL=struct('Type',Type,'Thickness',Thickness,'Bmax',Bmax,'miur',miur,'Kst',Kst,'Kh',Kh,...
                   'Kex',Kex,'Ke',Ke,'rhobi',rhobi,'Vst',Vst,'Wst',Wst,'Vro',Vro,'Wro',Wro);
%--------------------------------------------------------------------------
% WINDING:
%--------------------------------------------------------------------------
% Conductor Type:
Type='Copper';
% Copper resistivity at 20°C:
rhocu= 1.724138e-8;                % Ohms-m
% Copper Density:
%Dcu= 8960;                        % kg/m3
Dcu= 8933;                         % kg/m3
% Temperature Coefficient:
Beta= 0.004041;
% Copper resistivity at 75°C:
rhocu75= rhocu*(1+Beta*(75-20));
% Copper resistibity at Tx (°C):
rhocux=rhocu;                       % Ohms-m
% Conductor Packing Factor (Stator Winding Factor):
kcp = 0.879539;
kcp=0.5;
% Conductor Area:
Ac=[];                         % m2
% Diameter Conductor:
dcu=1.8;                       % mm
% Conductor Weigth:
Wcu=[];                        % kg
% Conductor volume:
Vcu=[];                        % m3
% End winding Geometry:
% Coil thickness:
bcu=[];                        % m
% Spacing between coils:
te=[];                         % m
alfa=[];                       % m
le3=[];                        % m
le2=0.06;                      % m
le1=0.02;                      % m
% winding type:
Winding=[];
winding_type=[];
% Winding Distribution:
WDT_OneLayer=[];
WDT_DoubleLayer=[];
DChain_OneLayer=[];
DChain_DoubleLayer=[];
% Other windings Wave, Double-Layer
PMSG.WINDING=struct('Type',Type,'Winding',Winding,'winding_type',winding_type,'WDT_OneLayer',WDT_OneLayer,...
    'DChain_OneLayer',DChain_OneLayer,'WDT_DoubleLayer',WDT_DoubleLayer,'DChain_DoubleLayer',DChain_DoubleLayer,...
    'rhocu',rhocu,'rhocu75',rhocu75,'rhocux',rhocux,'Dcu',Dcu,'Beta',Beta,'kcp',kcp,'Ac',Ac,'dcu',dcu,'Wcu',Wcu,...,
    'Vcu',Vcu,'bcu',bcu,...
    'alfa',alfa,'le1',le1,'le2',le2,'le3',le3);
%--------------------------------------------------------------------------
% SHOE:
%--------------------------------------------------------------------------
% Slot Opening:
%ws=0.0023;                     % m
ws=[];                          %
% Slot Fraction (as=ws/Ts):
as=0.25;
% Shoe Depth Fraction (d1+d2)/wsb:
asd= 0.5238;                   %
% Other windings Wave, Double-Layer
PMSG.SHOE=struct('ws',ws,'as',as,'asd',asd);
%--------------------------------------------------------------------------
% CONSTANTS:
%--------------------------------------------------------------------------
% Distribution Factor:
kd=[];
% Pole-Pitch Factor:
kp=[];
% Skew Factor:
ks=[];
% Conductor Packing Factor:
kcp=0.8;
PMSG.FACTORS=struct('kd',kd,'kp',kp,'ks',ks,'kcp',kcp);
%--------------------------------------------------------------------------
% AIR-GAP FLUX:
%--------------------------------------------------------------------------
% Air-gap Flux:
Phig=[];                    % Wb
% Back-iron flux:
Phibi=[];                   % Wb
% Air-gap Flux Density:
Bg=[];                      % T
% Air-gap Flux Density:
Bgx=[];                     % Tx
% Back-iron flux density:
Bbi=[];                     % T
PMSG.FLUX=struct('Phig',Phig,'Phibi',Phibi,'Bg',Bg,'Bgx',Bgx,'Bbi',Bbi);
%--------------------------------------------------------------------------
% ELECTRICAL PARAMETERS:
%--------------------------------------------------------------------------
% Number of Turns per Slot:
ns=[];                      % Number of Turns per slot
nsx=[];                      % Number of Turns per slot (Dependance of Tx)
% emax:
Emax=[];                    % Volts
Emaxx=[];                   % Volts (Dependance of Tx)
% Total Slot Current:
Is=[];                      % Amp
Isx=[];                     % Amp (in dependence with Tx)
% Conductor current:
Ic=[];                      % A
Icx=[];                     % A (in dependence with Tx)
% Phase Current:
Iph=[];                     % Amp
Iphx=[];                     % Amp (in dependence with Tx)
% Maximum current density:
Jmax=4.5;                   % A/mm2
% Slot and conductor current density:
Jc=[];                      % A/m2
Jcx=[];                      % A/m2 (in dependence with Tx)
% Bsmax:
Bsmax=[];                    % Wb/m2
Bsmaxx=[];                   % Wb/m2  (in dependence with Tx)
% Phase Resistance at temperature Tx:
Rph=[];                      % Ohms (20 °C)
Rphx=[];                     % Ohms
% Phase Resistance (75°C):
Rph75=[];                    % Ohms
% Slot Resistance;
Rs=[];                      % Ohms
Rsx=[];                     % Ohms  (in dependence with Tx)
% End-turn resistance:
Re=[];                      % Ohms
Rex=[];                     % Ohms  (in dependence with Tx)
% Phase Inductance:
Lph=[];                     % Henrios
% Phase Inductance:
Lphx=[];                    % Henrios  (in dependence with Tx)
% Air-Gap Inductance:
Lg=[];                      % Henrios
% Air-Gap Inductance:
Lgx=[];                     % Henrios  (in dependence with Tx)
% Slot-Leakage Inductance:
Ls=[];                      % Henrios
Lsx=[];                     % Henrios  (in dependence with Tx)
% End-Turn Inductance per Slot:
Le=[];                      % Henrios
Lex=[];                     % Henrios  (in dependence with Tx)
% Phase Reactance:
Xph=[];                     % Ohms
Xphx=[];                    % Ohms     (in dependence with Tx)
% Air-gap Reactance:
Xg=[];                      % Ohms
Xgx=[];                     % Ohms     (in dependence with Tx)
% Slot-Leakage Reactance:
Xs=[];                      % Ohms
Xsx=[];                     % Ohms     (in dependence with Tx)
% End-turn Reactance
Xe=[];                      % Ohms
Xex=[];                     % Ohms     (in dependence with Tx)
% STEADY STATE PARAMETERS IN dq0-Axis:
Lq=[];                      % H
Lqx=[];                     % H     (in dependence with Tx)
Ld=[];                      % H
Ldx=[];                     % H     (in dependence with Tx)
L0=[];                      % H
L0x=[];                     % H     (in dependence with Tx)
Xq=[];                      % Ohms
Xqx=[];                     % Ohms  (in dependence with Tx)
Xd=[];                      % Ohms
Xdx=[];                     % Ohms  (in dependence with Tx)
X0=[];                      % Ohms
X0x=[];                     % Ohms  (in dependence with Tx)
% Constants:
cte=[];                     % For evaluating Erms,
ctex=[];                    % To calculate Erms, Tx dependence
PMSG.ELECTRICAL=struct('ns',ns,'nsx',nsx,'Emax',Emax,'Emaxx',Emaxx,'Is',Is,'Isx',Isx,'Ic',Ic,'Icx',Icx,'Iph',Iph,'Iphx',Iphx,...
                       'Jmax', Jmax,'Jc',Jc,'Jcx',Jcx,'Bsmax',Bsmax,'Bsmaxx',Bsmaxx,'Rphx',Rphx,'Rph',Rph,'Rs',Rs,'Rsx',Rsx,'Re',Re,...
                       'Rex',Rex,'Lph',Lph,'Lphx',Lphx,'Ls',Ls,'Lsx',Lsx,'Lg',Lg,'Lgx',Lgx,'Le',Le,'Lex',Lex,'Xph',Xph,'Xphx',Xphx,'Xg',Xg,...
                       'Xgx',Xgx,'Xs',Xs,'Xsx',Xsx,...
                       'Xe',Xe,'Xex',Xex,'Lq',Lq,'Lqx',Lqx,'Ld',Ld,'Ldx',Ldx,'L0',L0,'L0x',L0x,'Xq',Xq,...
                       'Xqx',Xqx,...
                       'Xd',Xd,'Xdx',Xdx,'X0',X0,'X0x',X0x,'cte',cte,'ctex',ctex);
%--------------------------------------------------------------------------
% THERMAL CONSTANTS:
%--------------------------------------------------------------------------
% Heat convection coefficient:
% Cooling air density:
rhoc=1.1;                                     % kg/m3
% Heat capacity of the cooling air:
kthc=1010;                                    % J/(kg°K)
% Total volumetric cooling air flow:
vair=0.1;                                     % m/s
% Height cooling air channel:
hair=1e-10;                                   % m
% Frame length:
hframe=8.822e-3;                              % m
% Stator yoke back:
asy= 60;                                      % W/(°K-m2)
% Airgap:
ag= 40;                                       % W/(°K-m2)
% End-shields:
aes= 25;                                      % W/(°K-m2)
% End winding:
aew= 25;                                      % W/(°K-m2)
% Rotor yoke back:
ary= 25;                                      % W/(°K-m2)
% Frame:
afr= 7.1;                                     % W/(°K-m2)
% Shaft:
% ???????
% Bearings:
% ???????
% Thermal conductivities:
% Air (considering laminar flow) (0,025-0,03):
lamair=0.025;                                  % W/(°K-m)
% Iron:
lamfe=38;                                      % W/(°K-m)
% Coil:
lamcoil=1.8;                                   % W/(°K-m)
% Copper:
lamcu=400;                                     % W/(°K-m)
% Insulation:
lami=0.2;                                      % W/(°K-m)
% Iron frame: (By assuming that Tamb= 30 °C and Temperature rise (DT)=50 °C):
% See Table 7.1 of Chapter 7: Thermal Design
% Lipo, Introduction to AC Machine Design.
lamf=7.1;                                      % W/(°K-m)
% Thermal conductivity of magnets:
lamm=9;                                        % W/(°K-m)
% Magnet glue (Pegamento magnético):
lamglue=0.7;                                   % W/(°K-m)
% GRP magnet protection:
lamgrp=0.2;                                    % W/(°K-m)
% Thickness of magnet glue:
lmglue=0.1e-3;                                 % m
% Thickness of magnet reinforcement:
lgrp=0.5e-3;                                   % m
% Coil insulation thickness:
lci=1e-3;                                      % m
% Air temperature:
Tair=30;                                       % °C
% Copper Temperature:                          % °C
Tcu=Tair;                                      % °C
% Magnet Temperature:
Tm=Tair;                                       % °C
% Gap Temperature:
Tgap=Tair;                                     % °C
% Stator Temperature:
Tst=Tair;                                      % °C
PMSG.THERMAL=struct('rhoc',rhoc,'kthc',kthc,'vair',vair,'hair',hair,'hframe',hframe,'asy',asy,'ag',ag,'aes',aes,'aew',aew,...
                   'ary',ary,'afr',afr,'lamfe',lamfe,'lamcoil',lamcoil,'lamair',lamair,'lamf',lamf,...
                   'lamcu',lamcu,'lami',lami,'lamm',lamm,'lamglue',lamglue,'lamgrp',lamgrp,'lmglue',lmglue,...
                   'lgrp',lgrp,'lci',lci,'Tair',Tair,'Tcu',Tcu,'Tm',Tm,'Tgap',Tgap,'Tst',Tst);
%--------------------------------------------------------------------------
% PERFORMANCE:
%--------------------------------------------------------------------------
% Eficiency:
eta=[];                     % 
% Output Power:
Pout=[];                   % Watts
% Input Power:
Pin=[];                    % Watts
% Ohmic Loss Power:
Pr=[];                     % W
% Core Loss:
Pcl=[];                    % W
% Heat-Density in W/m2 leaving the slot conductor:
qr=[];                      % W/m2
% Maximum Heat-Density to be removed:
qrt=[];

PMSG.PERFORMANCE=struct('eta',eta,'Pout',Pout,'Pin',Pin,'Pr',Pr,'Pcl',Pcl,'qr',qr,'qrt',qrt);
%--------------------------------------------------------------------------
% SAVE PARAMETERS PMSG:
%--------------------------------------------------------------------------
% savefile='PMSG.mat';
% save(savefile,'PMSG');
