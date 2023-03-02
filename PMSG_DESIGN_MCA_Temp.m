% PMSG DESIGN PROGRAM0.0515

% The basis of this program is the Magnetic Circuit Analysis (MCA)

% Date: July 8th, 2017
clc;
close all;
clear;
%Tm=am*Tp;
%--------------------------------------------------------------------------
%PMSG PARAMETERS:
%--------------------------------------------------------------------------
[PMSG]=PMSG_PARAMETERS();
%load PMSG.mat
%--------------------------------------------------------------------------
% Topological Constraints:
%--------------------------------------------------------------------------
[PMSG]=TOPOLOGICAL_CONSTRAINT(PMSG);
%--------------------------------------------------------------------------
% PMSG Geometry:
%--------------------------------------------------------------------------
[PMSG]=PMSG_GEOMETRY(PMSG);
%PMSG.GEOMETRY.am=0.6711;            % %
%--------------------------------------------------------------------------
% Magnet Sizing:
%--------------------------------------------------------------------------
[PMSG]=PMSG_MAGNET(PMSG);
T=[50 30 70 90 110];                % °C
%Tm(1)=T(1);
%Tmm=50;      % °C
for i=1:5
%--------------------------------------------------------------------------
% Flux Calculations:
%--------------------------------------------------------------------------
% Fake:
%Tm(i)=T(i);
PMSG.THERMAL.Tm=T(i);
[PMSG]=PMSG_FLUX(PMSG);
Br(i)=PMSG.MAGNET.Br;
Brx(i)=PMSG.MAGNET.Brx;
Bg(i)=PMSG.FLUX.Bg;
Bgx(i)=PMSG.FLUX.Bgx;
Tm(i)=PMSG.THERMAL.Tm;
%--------------------------------------------------------------------------
% FACTORS:
%--------------------------------------------------------------------------
[PMSG]=PMSG_FACTORS(PMSG);
%--------------------------------------------------------------------------
% ELECTRICAL PARAMETER:
%--------------------------------------------------------------------------
[PMSG]=PMSG_ELECTRICAL(PMSG);
Rph(i)=PMSG.ELECTRICAL.Rph;
Rphx(i)=PMSG.ELECTRICAL.Rphx;
Xph(i)=PMSG.ELECTRICAL.Xph;
Xphx(i)=PMSG.ELECTRICAL.Xphx;
Tcu(i)=PMSG.THERMAL.Tcu;
rho(i)=PMSG.WINDING.rhocu;
rhox(i)=PMSG.WINDING.rhocux;
Emax(i)=PMSG.ELECTRICAL.Emax;
Emaxx(i)=PMSG.ELECTRICAL.Emax;
ns(i)=PMSG.ELECTRICAL.ns;
nsx(i)=PMSG.ELECTRICAL.nsx;
np(i)=PMSG.ELECTRICAL.np;
npx(i)=PMSG.ELECTRICAL.npx;
Iph(i)=PMSG.ELECTRICAL.Iph;
Iphx(i)=PMSG.ELECTRICAL.Iphx;
%--------------------------------------------------------------------------
% WINDING PARAMETER:
%--------------------------------------------------------------------------
[PMSG]=PMSG_WINDING(PMSG);
%--------------------------------------------------------------------------
% THERMAL PARAMETER:
%--------------------------------------------------------------------------
[PMSG]=PMSG_THERMAL(PMSG);
end
%--------------------------------------------------------------------------
% PMSG PERFORMANCE:
%--------------------------------------------------------------------------
[PMSG]=PMSG_PERFORMANCE(PMSG);
%--------------------------------------------------------------------------
% PLOTS: Temperatures:
%--------------------------------------------------------------------------
PMSG_PLOTS(PMSG);
%--------------------------------------------------------------------------
% PLOTS: Broad Range Performance:
%--------------------------------------------------------------------------
% Rotor speeds:
[PMSG]=PMSG_TEST(PMSG);
%--------------------------------------------------------------------------
% SAVE PARAMETERS PMSG:
%--------------------------------------------------------------------------
%disp(PMSG)
savefile='PMSG1.mat';
save(savefile,'PMSG');
