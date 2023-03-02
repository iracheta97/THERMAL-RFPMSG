function PMSG=PMSG_TEST(PMSG)
% Help: This function do the broad range performance test of the RF_PMSG
% Created by: Dr. Reynaldo Iracheta Cortez
% Date: June 11th, 2021
% Inputs:
% PMSG: Structure Data of the PMSG
% Outputs:
% PERFORMANCE TEST PLOTS
% -------------------------------------------------------------------------
% Input Data:
%--------------------------------------------------------------------------
% Number of phases:
Nph=PMSG.CONSTRAINTS.Nph;
% Resistance:
Rph=PMSG.ELECTRICAL.Rph;              % Ohms
Rphx=PMSG.ELECTRICAL.Rphx;            % Ohms
% Inductance:
Lph=PMSG.ELECTRICAL.Lph;              % H
Lphx=PMSG.ELECTRICAL.Lphx;            % H
% Reactance (60 Hz):
Xph=PMSG.ELECTRICAL.Xph;              % Ohms
Xphx=PMSG.ELECTRICAL.Xphx;            % Ohms
% Number of series turns:
ns=PMSG.ELECTRICAL.ns;
nsx=PMSG.ELECTRICAL.nsx;
% Number of parallel turns:
np=PMSG.ELECTRICAL.np;
npx=PMSG.ELECTRICAL.npx;
% Constant:
cte=PMSG.ELECTRICAL.cte;
ctex=PMSG.ELECTRICAL.ctex;
% Nm: Number of magnets
Nm= PMSG.CONSTRAINTS.Nm;
% Voltage vs. Power at f=60 Hz;
fe=PMSG.fe.Data;                     % Hz
nrpm=120*fe/Nm;                      % rpm
wm60=2*pi*nrpm/60;                   % rad/s
% Erms60=cte*wm60*ns/sqrt(2);
% Erms60x=ctex*wm60*nsx/sqrt(2);
Erms60=cte*wm60*ns;
Erms60x=ctex*wm60*nsx;
% Vector of current in %:
P3=PMSG.Pout.Data;                   % W
Vn=PMSG.Emax.Data;                   % V
In= P3/(sqrt(3)*Vn);                 % A
% Vector of Resistances:
Rph=PMSG.ELECTRICAL.Rph;             % Ohms
Rphx=PMSG.ELECTRICAL.Rphx;           % Ohms
% Steel:
Bmax=PMSG.STEEL.Bmax;
kh=PMSG.STEEL.Kh;
ke=PMSG.STEEL.Ke;
Vst=PMSG.STEEL.Vst;                 % m3
%--------------------------------------------------------------------------
% Parametric Study:
%--------------------------------------------------------------------------
% Operation speed range n(rpm):
n=[0 25 50 100 150 200 225 250 300 325 350 375 400 425 450 475];     % rpm
% Operation load range (1-100 %);
LP=linspace(1,100,length(n));
% To find data:
m1=find(n==225);
m2=find(n==max(n));
LP50=LP(m1);
LP100=LP(m2);
n225=n(m1);
n475=n(m2);
% abs(P):
AP=LP*P3/100;                            % W
% Load Resistance:
RR=(Vn*Vn)./(AP);
% Rank of resistances R(Ohms):
RR1x=RR+Rphx+j*Xphx;
RR1=RR+Rph+j*Xph;
% Rank of currents I(Amp):
Ix=Erms60x/sqrt(3)./RR1x;                 % A (In dependence with Tx)
I=Erms60/sqrt(3)./RR1;                    % A 
% Power:
Px=3*Ix.*Ix.*RR1x;                       % W (In dependence with Tx)
P=3*I.*I.*RR1;                           % W
% Voltage profile for different percentage loads:
VTx=sqrt(3)*(abs(Erms60x)/sqrt(3)-abs(Ix*(Rphx+j*Xphx))); % W (In dependence with Tx)
VT=sqrt(3)*(abs(Erms60)/sqrt(3)-abs(I*(Rph+j*Xph)));
% Number of samples:
N=length(n);
for i=1:N
    % wm (rad/s):
    wm(i)=2*pi*n(i)/60;
    % fee:
    fee(i)=(Nm/120)*n(i);
    % we (rad/s):
    we(i)= (Nm/2)*wm(i);
    %we1(i)=2*pi*fee(i);
    % Internal Voltage: Erms (V):
    %Erms1(i)=cte*wm(i)*ns/sqrt(2)
    Erms1(i)=cte*wm(i)*ns;
    Erms1x(i)=ctex*wm(i)*ns;
    % Current at different speeds (A):
    % RLoad= 14.52;                         % Ohms;
    RLoad= 6.25;                            % Ohms
    % RLoad= Erms1(i)*Erms1(i)/P3;          % Ohms
    % Rated Current:
    ILx(i)= Erms1(i)/sqrt(3)/(RLoad+Rphx+j*2*pi*fee(i)*Lph);           % 100%
    IL(i)= Erms1(i)/sqrt(3)/(RLoad+Rph+j*2*pi*fee(i)*Lph);             % 100%
    % Voltage profile at Full-load:
    VT1x(i)=abs(Erms1x(i)-abs((ILx(i)*(Rphx+j*2*pi*fee(i)*Lphx))));
    VT1(i)=abs(Erms1(i)-abs((IL(i)*(Rph+j*2*pi*fee(i)*Lph))));
    % Probably the calculation of phase inductance is wrong,
    %VT1(i)=0*Erms1(i)-abs((In*(Rph+j*Xph)));
    % VT2 (V) for broad ranges:
    for p=1:length(LP)
        % Rated Current 2:
        IL2x(p)=(Erms1x(i)/sqrt(3))/(RR(p)+Rphx+j*2*pi*fee(i)*Lphx);               % 100%
        IL2(p)=(Erms1(i)/sqrt(3))/(RR(p)+Rph+j*2*pi*fee(i)*Lph);                 % 100%
        % VLOAD:
        VT2x(i,p)=(Erms1x(i)-abs(IL2x(p)*(Rphx+j*2*pi*fee(i)*Lphx)));            % Volts
        VT2(i,p)=(Erms1(i)-abs(IL2(p)*(Rph+j*2*pi*fee(i)*Lph)));              % Volt
        % Copper Losses: W
        PCUx(i,p)=Nph*abs(IL2x(p).*IL2x(p)*Rphx);                                % Watts
        PCU(i,p)=Nph*abs(IL2(p).*IL2(p)*Rph);                                    % Watts
        % Core Losses: W
        PCLx(i,p)=(kh*fee(i)*Bmax^2+ke*fee(i).^2*Bmax^2+ke*fee(i).^1.5*Bmax^1.5)*Vst;          % Watts
        PCL(i,p)=(kh*fee(i)*Bmax^2+ke*fee(i).^2*Bmax^2+ke*fee(i).^1.5*Bmax^1.5)*Vst;           % Watts
        %PCL(i,p)=(kh*60*Bmax^2+ke*60.^2*Bmax^2+ke*60.^1.5*Bmax^1.5)*Vst;
        % Efficiency:
        Poutx(i,p)=sqrt(3)*VT2x(i,p)*IL2x(p);                              % Watts
        Pout(i,p)=sqrt(3)*VT2(i,p)*IL2(p);                              % Watts
        % Input power (W):
        Pinx(i,p)=Poutx(i,p)+PCUx(i,p)+PCLx(i,p)+0.01*Poutx(i,p);          % Watts
        Pin(i,p)=Pout(i,p)+PCU(i,p)+PCL(i,p)+0.01*Pout(i,p);               % Watts
        % Efficiency:
        etax(i,p)=100*abs(Poutx(i,p))/abs(Pinx(i,p));                      % %
        eta(i,p)=100*abs(Pout(i,p))/abs(Pin(i,p));                         % %
        % Voltage Regulation (%):
        Regx(i,p)=100*(VT2x(i,1)-VT2x(i,p))/VT2x(i,p);                     % %
        Reg(i,p)=100*(VT2(i,1)-VT2(i,p))/VT2(i,p);                         % %
    end
end
% Save Data:
PMSG.TEST.Erms1=Erms1;             % V
PMSG.TEST.Erms1x=Erms1x;           % V
% Voltage across terminals:
PMSG.TEST.VT1=VT1;                 % V
PMSG.TEST.VTx=VT1x;                % V (Tx)
% Load Test:
Erms=ones(1,length(I))*Erms60;
Ermsx=ones(1,length(I))*Erms60x;   % In dependence with Tx
PMSG.TEST.Erms=Erms;               % V
PMSG.TEST.Ermsx=Ermsx;             % V
PMSG.TEST.VTL=VT;                  % V
PMSG.TEST.VTLx=VTx;                % V
% 3D Voltage profile:
PMSG.TEST.VT3D=VT2;                % V
PMSG.TEST.VT3Dx=VT2x;              % V
% Absolute error:
EV3D=VT2-VT2x;
PMSG.TEST.EV3D=VT2-VT2x;            % V
% Copper losses:
PMSG.TEST.Pcu=PCU;                 % W
PMSG.TEST.Pcux=PCUx;               % W
% Absolute error:
EP3D=PCUx-PCU;
PMSG.TEST.EP3D=EP3D;               % W
% Core losses:
PMSG.TEST.Pcl=PCL;                 % W
PMSG.TEST.Pclx=PCLx;               % W
% Absolute error:
EPCL=PCLx-PCL;
PMSG.TEST.EPCL=EPCL;               % W
%--------------------------------------------------------------------------
% PLOT: Voltage profile
%--------------------------------------------------------------------------
figure(1)
plot(n,Erms1, 'LineStyle','-','LineWidth',3,'color',[1 0 0]);
hold on;
plot(n,abs(Erms1x),'LineStyle','-.','LineWidth',3,'color',[0 0 0]);
hold off;hold on;
plot(n,abs(VT1),'LineStyle','-','LineWidth',3,'color',[0 0 1]);
hold off; hold on;
plot(n,abs(VT1x),'LineStyle','--','LineWidth',3,'color',[0 0 0]);
hold off;
%xlim([0 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('V_L_L (V)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
% NL = sprintf('No load, T= %.1f °C',PMSG.THERMAL.Tair);
% NLx = sprintf('No load, T= %.1f °C',PMSG.THERMAL.Tcu);
% FL= sprintf('Full-Load, T= %.1f °C',PMSG.THERMAL.Tair);
% FLx = sprintf('Full load, T= %.1f °C',PMSG.THERMAL.Tcu);
NL = sprintf('Vacío, T_c_u= %.1f °C',PMSG.THERMAL.Tair);
NLx = sprintf('Vacío, T_c_u= %.1f °C',PMSG.THERMAL.Tcu);
FL= sprintf('Plena carga, T_c_u= %.1f °C',PMSG.THERMAL.Tair);
FLx = sprintf('Plena carga, T_c_u= %.1f °C',PMSG.THERMAL.Tcu);
h=legend(NL,NLx,FL,FLx,'Location','Northwest');
%h.Font Size = 15;
lgd.FontWeight = 'normal';
%set(h,'LineStyle','-','LineWidth',2);
grid
%--------------------------------------------------------------------------
figure(2)
%--------------------------------------------------------------------------
% v2= [202.301 217.90 224.61 227.26 234.86 240.47 243.91 246.21 247.82 249.03 249.94];
% x2=100*[12.20 10.88 10.02 9.62 8.22 6.89 5.91 5.16 4.57 4.1 3.72]/10;
plot(LP,(Erms), 'LineStyle','-','LineWidth',3,'color',[1 0 0]);
hold on;
plot(LP,abs(Ermsx),'LineStyle','-.','LineWidth',3,'color',[0 0 0]);
hold off;hold on;
plot(LP,(VT), 'LineStyle','-','LineWidth',3,'color',[0 0 1]);
hold off; hold on;
plot(LP,abs(VTx),'LineStyle','--','LineWidth',3,'color',[0 0 0]);
hold off;
% hold on;
% plot(x2,v2, 'LineStyle','-','LineWidth',3,'color',[0 0 1]);
% hold off;
xlim([0 100]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('V_L_L (V)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
% xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
% 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
% NL = sprintf('No load, T= %.1f °C',PMSG.THERMAL.Tair);
% NLx = sprintf('No load, T= %.1f °C',PMSG.THERMAL.Tcu);
% FL= sprintf('Full-Load, T= %.1f °C',PMSG.THERMAL.Tair);
% FLx = sprintf('Full load, T= %.1f °C',PMSG.THERMAL.Tcu);
NL = sprintf('Vacío, T_c_u= %.1f °C',PMSG.THERMAL.Tair);
NLx = sprintf('Vacío, T_c_u= %.1f °C',PMSG.THERMAL.Tcu);
FL= sprintf('Plena carga, T_c_u= %.1f °C',PMSG.THERMAL.Tair);
FLx = sprintf('Plena carga, T_c_u= %.1f °C',PMSG.THERMAL.Tcu);
h=legend(NL,NLx,FL,FLx,'Location','Southwest');
h.FontSize = 15;
lgd.FontWeight = 'normal';
%set(h,'LineStyle','-','LineWidth',2);
grid
%--------------------------------------------------------------------------
% 3D Voltage profile:
%--------------------------------------------------------------------------
figure(3)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
%zlim([10,90]); 
colormap(flipud(colormap))
surf((LP), n, (VT2x),'FaceLighting','gouraud','LineWidth',0.5);view(50,20);
%view(40,20);
y.EdgeColor='k';
y.LineStyle='-';
y.LineWidth=1;
y.MarkerEdgeColor='k';
y.MarkerFaceColor=[.49 1 .63];
y.MarkerSize=1;
colorbar;
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
VV = sprintf('V_L_L',PMSG.THERMAL.Tcu);
zlabel(VV,'fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',-20,...
'HorizontalAlignment','center');
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',15,...
'HorizontalAlignment','center');
E1 = sprintf('V_L_L= %.1f V,\nS_r= %.1f rpm',[VT2x(m1,1) n(m1)]);
E2 = sprintf('V_L_L= %.1f V,\nS_r= %.1f rpm',[VT2x(m1,end) n(m1)]);
E3 = sprintf('V_L_L= %.1f V,\nS_r= %.1f rpm',[VT2x(m2,end) n(m2)]);
E4 = sprintf('V_L_L= %.1f V,\nS_r= %.1f rpm',[VT2x(m2,1) n(m2)]);
% text(-7.1418808247231, 84.1813212832822, 419.181124525222,E1,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(70.3683341329493, 367.614346135651, 73.1196523848566,E2,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(60.470916558681, 493.62821286088, 478.835294712671,E3,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(1.10442629664976, 277.423390751863, 661.023028497584,E4,'fontname','normal','FontSize',12,'BackgroundColor','white');
x1=[0.276091081593928,0.324003795066414];
y1=[0.602436323366556,0.573643410852708];
x2=[0.425047438330171,0.496204933586338];
y2=[0.823920265780731,0.85492801771871];
x3=[0.769449715370019,0.810246679316888];
y3=[0.673311184939092,0.546327057954965];
x4=[0.666982922201138,0.636148007590133];
y4=[0.321889996308601,0.366925064599481];
annotation('textarrow',x1,y1,'String',E1,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x2,y2,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x3,y3,'String',E3,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x4,y4,'String',E2,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
grid;grid;
%Call this to rotate the labels 
% rotate_labels(gca);
% %Optional Step. Cause labels to chage as figure is rotated 
% call_when_rotates = @(~,~,hax)(rotate_labels(hax)); 
% hrot = rotate3d; 
% set(hrot,'ActionPreCallBack',... 
%     @(hfig,event_data)(set(hfig,'WindowButtonMotionFcn',{call_when_rotates event_data.Axes}))); 
% set(hrot,'ActionPostCallBack',... 
%     @(hfig,~)(set(hfig,'WindowButtonMotionFcn','')));
%--------------------------------------------------------------------------
% ABSOLUTE ERROR 3D: VT(20-°C)-VT(80-°C)
%--------------------------------------------------------------------------
figure(4)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
surf((LP), n, EV3D,'FaceLighting','gouraud','LineWidth',0.5);view(50,20);
y.EdgeColor='k';
y.LineStyle='-';
y.LineWidth=1;
y.MarkerEdgeColor='k';
y.MarkerFaceColor=[.49 1 .63];
y.MarkerSize=1;
colorbar;
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
zlabel('\Delta V_L_L','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',-20,...
'HorizontalAlignment','center');
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',15,...
'HorizontalAlignment','center');
E1 = sprintf('\\DeltaV_L_L= %.1f V,\nS_r= %.1f rpm',[EV3D(m1,1) n(m1)]);
E2 = sprintf('\\DeltaV_L_L= %.1f V,\nS_r= %.1f rpm',[EV3D(m1,end) n(m1)]);
E3 = sprintf('\\DeltaV_L_L= %.1f V,\nS_r= %.1f rpm',[EV3D(m2,end) n(m2)]);
E4 = sprintf('\\DeltaV_L_L= %.1f V,\nS_r= %.1f rpm',[EV3D(m2,1) n(m2)]);
% text(-7.1418808247231, 84.1813212832822, 419.181124525222,E1,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(70.3683341329493, 367.614346135651, 73.1196523848566,E2,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(60.470916558681, 493.62821286088, 478.835294712671,E3,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(1.10442629664976, 277.423390751863, 661.023028497584,E4,'fontname','normal','FontSize',12,'BackgroundColor','white');
x1=[0.276091081593928,0.324003795066414];
y1=[0.602436323366556,0.573643410852708];
x2=[0.425047438330171,0.496204933586338];
y2=[0.823920265780731,0.85492801771871];
x3=[0.769449715370019,0.810246679316888];
y3=[0.673311184939092,0.546327057954965];
x4=[0.666982922201138,0.636148007590133];
y4=[0.321889996308601,0.366925064599481];
annotation('textarrow',x1,y1,'String',E1,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x2,y2,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x3,y3,'String',E3,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x4,y4,'String',E2,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
grid;grid;
%--------------------------------------------------------------------------
% Voltage drop contour map:
%--------------------------------------------------------------------------
% figure(5)
% colormap GRAY;
% % colormap hot;
% DVE=100*EV3D./abs(VT2);
% contourf(LP,n,DVE,[0  1 3 5]);
% xlim([0 100]);
% ylim([25 400]);
% set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
% ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
%  'HorizontalAlignment','center');
% % ylim([-150 200]);
% xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
% 'HorizontalAlignment','center');
% text(60,350,'\bf0.3% <\DeltaP_c_u \leq 0.5%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
% text(29,200,'\bf0.1% <\DeltaP_c_u \leq 0.3%\rm','FontSize',15,'BackgroundColor','white');
% %text(12,150,'\bf0.1% <\DeltaP_c_u \leq 0.2%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
% text(1,55,'\bf\DeltaP_c_u \leq 0.1%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
% %text(78,372,'\bf\eta > 95%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%--------------------------------------------------------------------------
% Power Losses:
%--------------------------------------------------------------------------
figure(6)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
xlim([0 500]);
surf(LP, n, (PCUx),'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
y.EdgeColor='k';
y.LineStyle='-';
y.LineWidth=1;
y.MarkerEdgeColor='k';
y.MarkerFaceColor=[.49 1 .63];
y.MarkerSize=1;
colorbar;
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
zlabel(' P_c_u (W)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',15,...
'HorizontalAlignment','center');
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',-20,...
'HorizontalAlignment','center');
%E1 = sprintf('P_c_u= %.1f W,\nS_r= %.1f rpm',[PCUx(m1,1) n(m1)]);
%E2 = sprintf('P_c_u= %.1f W,\nS_r= %.1f rpm',[PCUx(m1,end) n(m1)]);
E3 = sprintf('P_c_u= %.1f W,\nS_r= %.1f rpm',[PCUx(m1,end) n(m1)]);
E4 = sprintf('P_c_u= %.1f W,\nS_r= %.1f rpm',[PCUx(m2,end) n(m2)]);
% text(-7.1418808247231, 84.1813212832822, 419.181124525222,E1,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(70.3683341329493, 367.614346135651, 73.1196523848566,E2,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(60.470916558681, 493.62821286088, 478.835294712671,E3,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(1.10442629664976, 277.423390751863, 661.023028497584,E4,'fontname','normal','FontSize',12,'BackgroundColor','white');
x1=[0.276091081593928,0.324003795066414];
y1=[0.602436323366556,0.573643410852708];
x2=[0.425047438330171,0.496204933586338];
y2=[0.823920265780731,0.85492801771871];
x3=[0.663021582733813 0.639601245000205];
y3=[0.54728370221328 0.500526674257429];
x4=[0.528345323741007 0.454676258992806];
y4=[0.833802816901409 0.920724346076459];

%annotation('textarrow',x1,y1,'String',E1,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
%annotation('textarrow',x2,y2,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x3,y3,'String',E3,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x4,y4,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
grid;grid;
%--------------------------------------------------------------------------
% Absolute error in Power Losses:
%--------------------------------------------------------------------------
figure(7)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
xlim([0 500]);
surf(LP, n, (EP3D),'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
y.EdgeColor='k';
y.LineStyle='-';
y.LineWidth=1;
y.MarkerEdgeColor='k';
y.MarkerFaceColor=[.49 1 .63];
y.MarkerSize=1;
colorbar;
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
zlabel('\DeltaP_c_u (W)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',15,...
'HorizontalAlignment','center');
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',-20,...
'HorizontalAlignment','center');
%E1 = sprintf('P_c_u= %.1f W,\nS_r= %.1f rpm',[PCUx(m1,1) n(m1)]);
%E2 = sprintf('P_c_u= %.1f W,\nS_r= %.1f rpm',[PCUx(m1,end) n(m1)]);
E3 = sprintf('\\DeltaP_c_u= %.1f W,\nS_r= %.1f rpm',[EP3D(m1,end) n(m1)]);
E4 = sprintf('\\DeltaP_c_u= %.1f W,\nS_r= %.1f rpm',[EP3D(m2,end) n(m2)]);
% text(-7.1418808247231, 84.1813212832822, 419.181124525222,E1,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(70.3683341329493, 367.614346135651, 73.1196523848566,E2,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(60.470916558681, 493.62821286088, 478.835294712671,E3,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(1.10442629664976, 277.423390751863, 661.023028497584,E4,'fontname','normal','FontSize',12,'BackgroundColor','white');
x1=[0.276091081593928,0.324003795066414];
y1=[0.602436323366556,0.573643410852708];
x2=[0.425047438330171,0.496204933586338];
y2=[0.823920265780731,0.85492801771871];
x3=[0.663021582733813 0.639601245000205];
y3=[0.54728370221328 0.500526674257429];
x4=[0.528345323741007 0.454676258992806];
y4=[0.833802816901409 0.920724346076459];

%annotation('textarrow',x1,y1,'String',E1,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
%annotation('textarrow',x2,y2,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x3,y3,'String',E3,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x4,y4,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
grid;grid;
%--------------------------------------------------------------------------
% Power losses contour map:
%--------------------------------------------------------------------------
figure(8)
colormap GRAY;
% colormap hot;
DPE=100*EP3D./abs(Poutx);
contourf(LP,n,DPE,[0  0.1 0.3 0.5]);
xlim([0 100]);
ylim([25 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
text(60,350,'\bf0.3% <\DeltaP_c_u \leq 0.5%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(29,200,'\bf0.1% <\DeltaP_c_u \leq 0.3%\rm','FontSize',15,'BackgroundColor','white');
%text(12,150,'\bf0.1% <\DeltaP_c_u \leq 0.2%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(1,55,'\bf\DeltaP_c_u \leq 0.1%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%text(78,372,'\bf\eta > 95%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%--------------------------------------------------------------------------
% Increment of Power losses contour map:
%--------------------------------------------------------------------------
figure(9)
colormap GRAY;
% colormap hot;
plossx=100*PCUx./abs(Poutx);
contourf(LP,n,plossx,[0 1 3 5]);
xlim([0 100]);
ylim([25 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
text(67,350,'\bf3% <P_c_u \leq 5%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(32,200,'\bf1% <P_c_u \leq 3%\rm','FontSize',15,'BackgroundColor','white');
%text(12,150,'\bf1% <P_c_u \leq 2%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(1,55,'\bfP_c_u \leq 1%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%text(78,372,'\bf\eta > 95%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');

figure(10)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
xlim([0 500]);
surf(LP, n, (PCLx),'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
%view(40,20);
y.EdgeColor='k';
y.LineStyle='-';
y.LineWidth=1;
y.MarkerEdgeColor='k';
y.MarkerFaceColor=[.49 1 .63];
y.MarkerSize=1;
colorbar;
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
VV = sprintf('P_C_L (W)',PMSG.THERMAL.Tcu);
zlabel(VV,'fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',15,...
'HorizontalAlignment','center');
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',-20,...
'HorizontalAlignment','center');
E1 = sprintf('P_C_L= %.1f W,\nS_r= %.1f rpm',[PCLx(m1,1) n(m1)]);
E2 = sprintf('P_C_L= %.1f W,\nS_r= %.1f rpm',[PCLx(m1,end) n(m1)]);
E3 = sprintf('P_C_L= %.1f W,\nS_r= %.1f rpm',[PCLx(m2,end) n(m2)]);
E4 = sprintf('P_C_L= %.1f W,\nS_r= %.1f rpm',[PCLx(m2,1) n(m2)]);
% text(-7.1418808247231, 84.1813212832822, 419.181124525222,E1,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(70.3683341329493, 367.614346135651, 73.1196523848566,E2,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(60.470916558681, 493.62821286088, 478.835294712671,E3,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(1.10442629664976, 277.423390751863, 661.023028497584,E4,'fontname','normal','FontSize',12,'BackgroundColor','white');
x1=[0.305818433060059, 0.32946508526933];
y1=[0.363333333333334, 0.366024363233662];
x2=[0.653336795637497, 0.631524279407946];
y2=[0.582334664891257, 0.518419884598314];
x3=[0.216100561647499, 0.141214228403316];
y3=[0.800685812258894, 0.702957565366481];
x4=[0.515645894624231, 0.440759561380048];
y4=[0.833261894556365, 0.836690955850836];

%annotation('textarrow',x1,y1,'String',E1,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x2,y2,'String',E2,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
%annotation('textarrow',x3,y3,'String',E3,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x4,y4,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
grid;grid;
%--------------------------------------------------------------------------
% Absolute error: EPCL
%--------------------------------------------------------------------------
figure(11)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
xlim([0 500]);
surf(LP, n, (EPCL),'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
%view(40,20);
y.EdgeColor='k';
y.LineStyle='-';
y.LineWidth=1;
y.MarkerEdgeColor='k';
y.MarkerFaceColor=[.49 1 .63];
y.MarkerSize=1;
colorbar;
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
VV = sprintf('\\Delta P_C_L (W)');
zlabel(VV,'fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',15,...
'HorizontalAlignment','center');
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',-20,...
'HorizontalAlignment','center');
E1 = sprintf('\\DeltaP_C_L= %.1f W,\nS_r= %.1f rpm',[EPCL(m1,1) n(m1)]);
E2 = sprintf('\\DeltaP_C_L= %.1f W,\nS_r= %.1f rpm',[EPCL(m1,end) n(m1)]);
E3 = sprintf('\\DeltaP_C_L= %.1f W,\nS_r= %.1f rpm',[EPCL(m2,end) n(m2)]);
E4 = sprintf('\\DeltaP_C_L= %.1f W,\nS_r= %.1f rpm',[EPCL(m2,1) n(m2)]);
% text(-7.1418808247231, 84.1813212832822, 419.181124525222,E1,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(70.3683341329493, 367.614346135651, 73.1196523848566,E2,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(60.470916558681, 493.62821286088, 478.835294712671,E3,'fontname','normal','FontSize',12,'BackgroundColor','white');
% text(1.10442629664976, 277.423390751863, 661.023028497584,E4,'fontname','normal','FontSize',12,'BackgroundColor','white');
x1=[0.305818433060059, 0.32946508526933];
y1=[0.363333333333334, 0.366024363233662];
x2=[0.653336795637497, 0.631524279407946];
y2=[0.582334664891257, 0.518419884598314];
x3=[0.216100561647499, 0.141214228403316];
y3=[0.800685812258894, 0.702957565366481];
x4=[0.515645894624231, 0.440759561380048];
y4=[0.833261894556365, 0.836690955850836];

%annotation('textarrow',x1,y1,'String',E1,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x2,y2,'String',E2,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
%annotation('textarrow',x3,y3,'String',E3,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
annotation('textarrow',x4,y4,'String',E4,'Color',[0.00,0.00,0.00],'LineStyle','-','LineWidth', 2.0,'HeadStyle','deltoid','HeadLength',12,'HeadWidth',12,'fontname','normal','FontSize',12);
grid;grid;
% Efficiency:
figure(12)
colormap Hot;
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
xlim([0 500]);
zlim([20 100]);
axis tight
surf(LP, n, eta,'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
set(gca,'fontname','tahoma','fontsize',15,'fontweight','light','LineWidth',2);

% Mapa de Contorno de la eficiencia:
figure(13)
colormap GRAY;
% colormap hot;
contourf(LP,n,eta,[0 60 80 90 95 100]);
xlim([0 100]);
ylim([25 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
text(60,250,'\bf90% <\eta \leq 95%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(47,75,'\bf80% <\eta \leq 90%\rm','FontSize',15,'BackgroundColor','white');
text(10,80,'\bf60% <\eta \leq 80%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(1,40,'\bf\eta \leq 60%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
text(78,372,'\bf\eta > 95%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%contour(DP,(ETA),MGL,30);
%contour3(DP,(ETA),MGL,30);

% Mapa de Contorno de la Regulación de Voltaje:
figure(14)
colormap GRAY;
% colormap hot;
contourf(LP,n,Regx,[0 5 20 50 70 100]);
xlim([0 100]);
ylim([25 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
%text(60,255,'\bf5%< V_r_e_g \leq10%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(47,145,'\bf5% < V_r_e_g \leq 20%\rm','FontSize',15,'BackgroundColor','white');
%text(40,65,'\bf5% < V_r_e_g \leq 2%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(5,60,'\bfV_r_e_g \leq 5%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
text(76,370,'\bfV_r_e_g > 20%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%contour(DP,(ETA),MGL,30);
%contour3(DP,(ETA),MGL,30);