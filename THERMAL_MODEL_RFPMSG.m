% This program solves the thermal model of a INDUCTION MACHINE
% Created by: Dr. Reynaldo IRacheta Cortez
% Date: June 26th,2020

clc;
clear;

% Variable Sleep:
s=2:-0.005:-1;

for i=1:length(s)
% Input parameters:
VLL=2400;                 % rms
% Conexión: 1: Estrella/Wye, 2: Delta
CX=1;
if CX==1    % Estrella%Wye:
    VLN= VLL/sqrt(3);
else
    VLN= VLL;
end
Pout=186.5e3;             % kW (250 HP)
% Copper resistivity (Ohm-m):
rhocu=1.68e-8;                      % Ohm-m

% Slip:
%s=0.03125;                          % slip
% Circuit Parameters:
% Stator:
R1=0.323;                           % Ohms
R2=0.786;                           % Ohms
Z1=R1+j*2.27+0.086;                 % Ohms
Z2= 1i*3.69+R2*s(i)+R2/s(i);            % Ohms
Z2= 1i*3.69+R2/s(i);            % Ohms
Z3= 1559*(j*84.7)/(1559+j*84.7);    % Ohms
Zeq=Z1+Z2*Z3/(Z2+Z3);
% Currents:
I1=(VLN/Zeq);                       % Amp
I2=I1*Z3/(Z2+Z3);                   % Amp
% Currents (Absolute values):
I11(i)=abs(I1);                        % Amp
I22(i)=abs(I2);                        % Amp

% Power losses analysis:
% Stator:
Pcu1=3*R1*I11(i).^2;                   % Watts

% Total Number of turns
Nt=240;
Ns=Nt;
% Number of circuits:
C=1;
% Poles:
P=8;                               %
% Pitch:
p=0.8;
% Pole pitch from the slot radial midpoint:
Tp=10.324*(2.54/100);              % m
% Slot pitch measured from the radial midpoint of the slots:
Ts=0.688*(2.54/100);               % m
% Longitud del núcleo bruta incluyendo el duct: lgc=li+2lo
lgc=9.25*(2.54/100);               % m
% Coil extension:
le2=1.25*(2.54/100);               % m
% Spacing between coils:
te=0.0625*(2.54/100);              % m
% Coil thickness:
bc=0.220*(2.54/100);               % m
% Area of coil:
Rfc=0.97;                          % Rounding factor of conductors
Ac=(0.129)*(0.204)*(Rfc)*(2.54/100)^2; % m2
% Length of the diagonal portion of the end winding:
p=0.8;                             % pitch
le3=(p*Tp/2)*Ts/(sqrt(Ts^2-(bc+te)^2));% m
% Mean length of a coil:
lc=2*lgc+4*le2+4*le3;               % m
% Stator resistance at 75 °C:
rhocu=1.72e-8;                      % Ohm-m
rhocu=2.1e-8;                       % Ohm-m
Rs=Nt*rhocu*lc/(C*Ac);
% Skin effect: Copper bar 75°C with 60 Hz
a0inv=(0.3707*2.54/100);            % m
a0=inv(a0inv);                      % m
% Thickness of a stator bar is 0.129''
d=0.129*2.54/100;                   % m
a0d=a0*d;
% AC RESISTANCE:
Rac=Rs*(1+(4/45)*(a0d)^4);          % m2

% Rotor resistance:
% Actual length including ducts:
lb1=9.25*(2.54/100);                % m (En el ejemplo 7, tiene la variable li) li=lb1
% Effective stator and rotor length:
les=8.67*2.54/100;                  % m
% Bar width:
b=0.375*2.54/100;                   % m
% Bar Depth:
d=0.5625*2.54/100;                  % m
% End ring area:
tbe= 1*2.54/100;                    % m
dbe=0.75*2.54/100;                  % m
% Rotor bar extension:
%lbe=0.375*2.54/100;                 % m
lbe=0.625*2.54/100;                 % m
% Average rotor pole pitch:
Tp2=8.604*2.54/100;                 % m
% Rotor slots:
S2=97;          % In PMSG S2=Ns2
% Length of the rotor bar:
lrb= 2*lbe+lb1;                     % m
% Active length of the rotor bars including the slight 
% additional length due to skee is:
lb=(lrb+(2/3)*tbe)/cos(2*pi/120);  % m
% Resistance of one rotor bar is, therefore:
rb= lb*rhocu/(b*d);                 % Ohms
% Tooth pitch at the middle of the end ring can be obtained from Tau_p2 as:
Tr2=Tp2*P/S2;                       % m
% Resistance of the end ring portion over one slot pitch is:
re=Tr2*rhocu/(tbe*dbe);             % Ohms
% Effective bar resistance:
rbe=rb+re/(2*sin(pi*P/(2*S2))^2);   % Ohms
% The resistance per rotor mesh:
rr=2*rbe;                           % Ohms
% Winding factor:
k1=0.91;                            % Dimensionless
% Skew factor:
ks1=0.9963;
% Rotor resistance referred to the stator:
rrs=12*k1^2*Ns^2*rr/S2;
% Skew effect:
rrss=ks1^2*rrs;                     % Ohms
% Lew End-winding
lew1=1.859*2.54/100;                    % m
lew=2*le2+2*sqrt(lew1^2+(p*Tp/2)^2);    % m
% Power dissipate in the stator slots:
Pslot1=Pcu1*lgc/(lgc+lew);              % W
Qslot1=Pslot1;
% Power dissipated in the each end-winding side:
Pew1=(Pcu1/2)*lew/(lgc+lew);            % W
Qew1=Pew1;
% Iron losses in the stator core:
Pcore1=630;                             % W
Pteeth1=290;                            % W
% Total heat generated within the stator core:
Qcore1=Pcore1+0.35*(Pcore1+Pteeth1);
Qteeth1=Pteeth1+0.7*(Pcore1+Pteeth1);
Qteeth2=0.7*(Pcore1+Pteeth1);
% Since the rotor experiences slip frequency during normal operations,
% The rotor core loss is typically neglected.
Qcore2=0;
% Copper losses:                         % W
Pcu2=3*R2*I22(i)^2;                      % W
% End ring resistance in the rotor:
rew2=re/(2*sin(pi*P/(2*S2))^2);          % Ohms
% Losses at the end ring:
Qew2=Pcu2*rew2/rbe;                      % Ohms
% Extra losses due to stray losses (Pérdidas de dispersión):
% It is considered 0.4% of rated output power as stray load loss
% The total heat loss injected into the equivalent circuit by the
% rotor conductor is:
Qslot2=Pcu2*rb/rbe+(0.4/100)*Pout;      % 4600 W

%--------------------------------------------------------------------------
% Thermal Resistances:
%--------------------------------------------------------------------------
% Stator tooth: Radial flow:
dss=2.2*(2.54/100);                     % m
S1=120;
S2=97;
li=9.25*(2.54/100);                     % m
li=8.5*(2.54/100);                     % m
ls=10*(2.54/100);                       % m
tts=0.256*(2.54/100);                   % m
tbs=0.369*(2.54/100);                   % m
% M-36 grade of silicon steel:
lamdairad=28;                           % W/°K-m
Rtrad1=dss/(lamdairad*S1*li*(tts+tbs)/2);
% Heat flow in the axial direction of the stator tooth:
lamdaiax=0.6;                           % W/°K-m
Rtax1=li/(lamdaiax*S1*dss*(tts+tbs)/2);
% Stator core in the radial direction:
Dos=31.5*(2.54/100);                    % m
Dbs=28.48*(2.54/100);                   % m
%Rcrad1=(1/(32*pi*lamdairad*li))*(Dos^2-3*Dbs^2+((4*Dbs^4)/(Dos^2-Dbs^2))*log(Dos/Dbs));
Rcrad1=(1/(8*pi*lamdairad*li))*(1-2*Dbs^2/(Dos^2-Dbs^2)+(4*Dbs^4/(Dos^2-Dbs^2)^2)*log(Dos/Dbs));
% For the stator core in the axial direction:
Rcax1=5*li/(6*pi*lamdaiax*(Dos^2-Dbs^2));
% For the rotor tooth, in the radial direction:
ds2=0.67*(2.54/100);                             % m
ttr=0.392*(2.54/100);                            % m
tbr=0.348*(2.54/100);                            % m
Rtrad2=ds2/(li*lamdairad*S2*li*(ttr+tbr)/2);      % W/°K
% For the rotor tooth in the axial direction:
Rtax2=li/(lamdaiax*S2*ds2*(ttr+tbr)/2);         % W/°K
% Stator conductors:
ls=9.25*2.54/100;                               % m
Acu1=12*0.129*0.204*(2.54/100)^2;               % m2
% Acu1=0.02553;                                 % m2
lamdacu=360;                                    % W/°K-m
Rcu1=ls/(lamdacu*S1*Acu1);                      % W/°K
% For the end/winding portion:
% Rew1=lew1/(lamdacu*S1*Acu1);                  % W/°K
Rew1=(lew/ls)*Rcu1;                             % W/°K
% Rotor conductors:
lamdacu=360;                                    % W/°K-m
Acu1=0.375*0.5625*(2.54/100)^2;                 % m2
Rcu2=ls/(lamdacu*S2*Acu1);                      % W/°K
% For the end-winding portion:
Aew=dbe*tbe;                                    % m2
LAEW= (lbe+tbe)/(Acu1)+0.71*(2.54/100)/Aew;     % 1/m
Rew2=LAEW/(lamdacu*S2);                         % W/°K
% Convection heat transfer from the end portions of the windind
% and core to the circulating air now be calculated:
% From one end-winding side, from equation (7.123):
% Insulation thickness:
tins1=0.078*(2.54/100);                          % m
% Thermal conductivity for insulation using synthetic resins:
lamdains1=0.25;                                  % W/(°K-m)
% The heat transfer coefficient for air passing over the end winding is:
vair1=2;                                         % m/s
aew1=20*vair1^0.6;                                % W/(°K-m2)
% Outside perimeter of the coil is:
dss=2.2*(2.54/100);                            % m
dcs=1.51*(2.54/100);                           % m
d2=0.317*(2.54/100);                           % m
d3=0.859*(2.54/100);                           % m
d4=0.0795*(2.54/100);                          % m
d5=0.858*(2.54/100);                           % m
d6=0.0858*(2.54/100);                          % m
bos=0.374*(2.54/100);                          % m
pcoil1=2*(d4/2+d5+d6+bos);                     % m
% Length of the end-winding side:
lend1=2*(le2+le3);
% Thermal resistance for one end-winding side is:
Rewcv1=(tins1/lamdains1+1/aew1)/(S1*pcoil1*lend1);   % °K/W
% Thermal resistance for one rotor and ring:
% Insulation thickness for the end ring is zero:
% The thermal heat transfer coefficient for the end ring is:
vair2=3;                                         % m/s
aew2=20*vair2^0.6;                               % W/(°K-m2)
% Mean diameter of the end ring:
Dor=24*(2.54/100);                              % m
d0r=0.04*(2.54/100);                            % m
d1r=0.047*(2.54/100);                           % m
d34r=0.572*(2.54/100);                          % m
d3r=0.5625*(2.54/100);                          % m
d4r=d34r-d3r;                                   % m
bsr=0.375*(2.54/100);                           % m
Dring=Dor-2*(d0r+d1r)-d34r;                     % m
lring=pi*Dring;                                 % m
% Perimeter of the ring:
pring=2*(tbe+dbe);                              % m
% The length and perimeter of the straight portions to the end of rotor
% bars:
pend2=2*(d3r+bsr);                              % m
lew2=(lbe+tbe);                                 % m
% Thermal resistance for one rotor bar and end ring:
Rewcv2=(1/aew2)/(S2*pend2*lew2+pring*lring);    % °K/W
% Thermal resistance of stator and rotor ducts:
% For stator:
aduct=30.3;                                           % W/(°K-m2)
lduct=3/8*(2.54/100);                                 % m
pcoil1=2*(d3+d4+d5+d6+bos);                           % m
Nd=2;                                                 % m2
Rduct1=(tins1/lamdains1+1/aduct)/(S1*pcoil1*lduct); % °K/W
% For Rotor:
% Assumption: Rotor bars are not insulated
pbar2=2*(d3r+bsr);
Rduct2=(1/aduct)/(S2*pbar2*lduct);                   % °K/W
% Insulation: Stator:
dss=(d3+d4+d5+d6)+0.08*2.54/100;                      % m
Rins1a= tins1/(lamdains1*S1*(bos+dss)*li);            % °K/W
% Air:
tair=0.035*(2.54/100);                                % m
lamdair=0.03;                                         % °K/W
Rins1b= tair/(lamdair*S1*(bos+dss)*li);               % °K/W
% Total thermal resistance between stator iron and stator conductors:
Rins1=Rins1a+Rins1b;                                   % °K/W
% Thermal resistance across the air gap:
g=0.04*(2.54/100);                                     % m
Dis=24.08*(2.54/100);                                  % m
Rag=g/(lamdair*pi*Dis*li);                             % °K/W
% Thermal resistance of stator iron also rises 
% The heat generated withing the stator iron also rises the
% outer surface of the stator, at which point heat is radiated and 
% convected away
% Radiated portion,
ar=7.1;                                                 % W/m2-°K
Arad=pi*Dos*li;                                         % m
Rfrad=1/(ar*Arad);                                      % °K/W
% Thermal resistance for convection of heat from the stator surface
% to the air must be calculated.
v=4;            
ac=7.8*v^0.78;                                          % W/(m2-°K)
Aconv=Arad;
Rfcv=1/(ac*Aconv);                                      % °K/W
% For the conbination of radiated and convected heat:
Rframe=Rfrad*Rfcv/(Rfrad+Rfcv);                         % °K/W
% Heat injected into the air in the inlet side of the machine:
Pdiss=0.5*(Pcu1+Pcu2+Qcore1+Qteeth1+Qteeth2);            % W
% Cubic meter per second:
CMS=0.236;                          % m3/s
% The air flow in cubic feet per minute (CFM):
DT=8.879e-4*5774/CMS;               % °K or °C
% Temperature of the cooling air at the midpoint of the machine:
Tair1=30;                           % °C
Tairm=Tair1+DT;
% Nodal analysis:
Taird=30;                           % °C
Tairl=Taird;                        % °C
Q1=Taird/(Rduct1+Rcu1/4);
Q2=Qew1+Tairl/(Rew1/2+Rewcv1);
Q3=Qew2+Tairl/(Rew2/2+Rewcv2);
Q4=(Qslot2+Qteeth2)/2+Taird/(Rcu2/4+Rduct2);
Q5= (Qslot1+Qteeth1)/2;
Q6=0;
Q7=Qcore1/2+Tairm/(2*Rframe);
Qh=[Q1;Q2;Q3;Q4;Q5;Q6;Q7];                  % W
% Matrix of thermal conductances:
G=zeros(7,7);                               % W/°K
% Self thermal conductances:
G(1,1)=1/(Rcu1/4+Rduct1)+1/(Rcu1/4+Rew1/2)+1/(2*Rins1);     % W/°K
G(2,2)=1/(Rcu1/4+Rew1/2)+1/(Rewcv1+Rew1/2);                 % W/°K
G(3,3)=1/(Rew2/2+Rewcv2)+1/(Rcu2/4+Rew2/2);                 % W/°K
G(4,4)= 1/(Rcu2/4+Rew2/2)+1/(Rcu2/4+Rduct2)+1/(2*Rag);      % W/°K
G(5,5)=1/(2*Rag)+1/(2*Rtrad1)+1/(2*Rins1);                  % W/°K
G(6,6)=1/(2*Rtrad1)+1/(2*Rcrad1);                           % W/°K
G(7,7)=1/(2*Rframe)+1/(2*Rcrad1);                           % W/°K
% Mutual thermal conductances:
G(1,2)=-1/(Rcu1/4+Rew1/2); 
G(1,5)=-1/(2*Rins1);
G(2,1)=G(1,2);
G(3,4)=-1/(Rew2/2+Rcu2/4);
G(4,3)=G(3,4);
G(4,5)=-1/(2*Rag);
G(5,1)=-1/(2*Rins1);
G(5,4)=G(4,5);
G(5,6)=-1/(2*Rtrad1);
G(6,5)=G(5,6);
G(6,7)=-1/(2*Rcrad1);
G(7,6)=G(6,7);
% Vector of Nodal temperatures:
T(:,i)=inv(G)*Qh;                                                % °C
end
%--------------------------------------------------------------------------
% Corrientes:
%--------------------------------------------------------------------------
figure(1)
plot(s,I11,'color',[0 .5 .5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
hold on;
plot(s,I22,'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(3,:),'color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% plot(s,I11(4,:),'color',[0 .5 .5],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(5,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(6,:),'color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(7,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% % %plot(s,[VL22],'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% % %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
hold off;
ylim([0 500]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Corriente (A)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('s','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('I_1, Estator','I_2,Rotor','Location','northeast','Orientation','vertical','NumColumns',1);
%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
%title(h,['Distribución de Temperaturas en la Máquina de Inducción']);
grid;
%--------------------------------------------------------------------------
% Temperatures:
%--------------------------------------------------------------------------
figure(2)
fe=60;               % Hz
ws=(2/P)*2*pi*fe;    % rads/s
Tmech=3*I22.^2*R2.*(1-s)./s;
plot(s,T(1,:),'color',[0 .5 .5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
hold on;
plot(s,T(2,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
plot(s,T(3,:),'color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
plot(s,T(4,:),'color',[0 .5 .5],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
plot(s,T(5,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
plot(s,T(6,:),'color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
plot(s,T(7,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
%plot(s,Tmech,'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
hold off;
ylim([0 5000]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Temperatura (°K)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('s','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('Devanado, Estator Centro', 'Devanado, Estator Extremo','Devanado,Rotor Extremo','Entrehierro, Rotor','Entrehierro, Estator', 'Diente, Estator','Yugo, Estator','Location','east','Orientation','vertical','NumColumns',2);
%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
title(h,['Distribución de Temperaturas en la Máquina de Inducción']);
grid;
%--------------------------------------------------------------------------
% Mechanical Torque (N-m):
%-------------------------------------------------------------------------
fe=60;               % Hz
ws=(2/P)*2*pi*fe;    % rads/s
Pmech=3*I22.^2*R2.*(1-s)./s;
figure(3)
plot(s,Pmech/1000,'color',[0 .5 .5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%hold on;
%plot(s,Tmech,'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(3,:),'color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% plot(s,I11(4,:),'color',[0 .5 .5],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(5,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(6,:),'color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(7,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% % %plot(s,[VL22],'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% % %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
%hold off;
%ylim([0 500]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Potencia (kW)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('s','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
%h = legend('I_1, Estator','I_2,Rotor','Location','northeast','Orientation','vertical','NumColumns',1);
%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
%title(h,['Distribución de Temperaturas en la Máquina de Inducción']);
grid;
figure(4)
plot(s,Tmech/1000,'color',[0 .5 .5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%hold on;
%plot(s,Tmech,'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(3,:),'color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% plot(s,I11(4,:),'color',[0 .5 .5],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(5,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(6,:),'color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(s,I11(7,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% % %plot(s,[VL22],'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% % %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
%hold off;
%ylim([0 500]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Torque (N-m)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('s','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
%h = legend('I_1, Estator','I_2,Rotor','Location','northeast','Orientation','vertical','NumColumns',1);
%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
%title(h,['Distribución de Temperaturas en la Máquina de Inducción']);
grid;