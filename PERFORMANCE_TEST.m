% Performance PLOTS of the TF-PMSG
% Created by: Dr. Reynaldo Iracheta Cortez,
% Date: January 30th, 2019

clc;
clear;
%--------------------------------------------------------------------------
% LOAD RF-PMSG DATA:
%--------------------------------------------------------------------------
load PMSG.mat
% Electrical Data:
% Phases:
Nph=PMSG.CONSTRAINTS.Nph;
% Resistance:
%Rph=PMSG.ELECTRICAL.Rph;            % Ohms
% Inductance:
Lph=(PMSG.ELECTRICAL.Lph);            % Ohms
% Reactance (60 Hz):
Xph=PMSG.ELECTRICAL.Xph;            % Ohms
% Number of seriesturns:
ns=PMSG.ELECTRICAL.ns;
% Constant:
cte=PMSG.ELECTRICAL.cte;
Nm= PMSG.CONSTRAINTS.Nm;


Erms1=250;
fe=60;
nrpm=120*fe/Nm;               % rpm
wm=2*pi*nrpm/60;            % rad/s
%ns1=floor(Erms1/(cte*wm))+1
%Erms=cte*wm*ns1
%Emax=Erms*sqrt(2)
%ns=floor(Erms/(cte*wm)+1)
%Erms=cte*wm*ns


% Number of parallel turns:
np=PMSG.ELECTRICAL.np;

% Nm:
Nm= PMSG.CONSTRAINTS.Nm;
% Voltage vs. Power at f=60 Hz;
fe=60;                        % Hz
nrpm=120*fe/Nm;               % rpm
wm60=2*pi*nrpm/60;            % rad/s
%Erms60=cte*wm60*ns/sqrt(2);
Erms60=(cte*wm60*ns)*sqrt(3)
%Erms60=PMSG.ELECTRICAL.Erms




% Vector of current in %:
P3=10e3;                      % W
Vn=Erms60;                       % V
In= P3/(sqrt(3)*Vn);          % A
% Vector of Resistances:
Rph=PMSG.ELECTRICAL.Rph;           % Ohms
% Steel:
Bmax=PMSG.STEEL.Bmax;
kh=PMSG.STEEL.Kh;
ke=PMSG.STEEL.Ke;
Vst=PMSG.STEEL.Vst;           % m3
Wst=PMSG.STEEL.Wst;
Wm=PMSG.MAGNET.Wm;
Wcu=PMSG.WINDING.Wcu;
Wrot=PMSG.STEEL.Wro;
Pcl=PMSG.PERFORMANCE.Pcl;
Pr=PMSG.PERFORMANCE.Pr;
Ps=PMSG.PERFORMANCE.Ps;
kc=PMSG.STEEL.Kc;
%--------------------------------------------------------------------------
% Parametric Study:
%--------------------------------------------------------------------------
% Rank of speeds n(rpm):
%n=linspace(1,400,N);                                                % rpm
n=[0 25 50 100 150 200 225 250 300 325 350 375 400 425 450 475 500];     % rpm
% Load Ranks (1-100 %);
LP=linspace(1,100,length(n));
% abs(P):
AP=LP*P3/100;                           % W
% Resistance:
RR=(Vn*Vn)./(AP);
% Rank of resistances R(Ohms):
RR1=RR+Rph+j*Xph;
% Rank of currents I(Amp):
I=Erms60/sqrt(3)./RR1;                  % A
% Power:
P=3*I.*I.*RR1;                          % W
% Voltage profile at different loads:
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
   % Erms1(i)=cte*wm(i)*ns/sqrt(2);
    Erms1(i)=(cte*wm(i)*ns)*sqrt(3);
    
    
    
    
    % Current at different speeds (A):
    % RLoad= 14.52;                           % Ohms;
    RLoad= 6.25;                            % Ohms
    % RLoad= Erms1(i)*Erms1(i)/P3;          % Ohms
    % Rated Current:
    ILoad(i)= Erms1(i)/sqrt(3)/(RLoad+Rph+j*2*pi*fee(i)*Lph);         % 100%
    % Voltage profile at Full-load:
    VT1(i)=sqrt(3)*(abs(Erms1(i)/sqrt(3)-abs((ILoad(i)*(Rph+j*2*pi*fee(i)*Lph)))));
    % Probably the calculation of phase inductance is wrong,
    %VT1(i)=0*Erms1(i)-abs((In*(Rph+j*Xph)));
    % VT2 (V) for broad ranges:
    for p=1:length(LP)
        % Rated Current 2:
        ILoad2(p)=(Erms1(i)/sqrt(3))/(RR(p)+Rph+j*2*pi*fee(i)*Lph);                            % 100%
        % VLOAD:
        VT2(i,p)=sqrt(3)*(Erms1(i)/sqrt(3)-abs(ILoad2(p)*(Rph+j*2*pi*fee(i)*Lph)));            % Volts
        % Copper Losses: W
        PCU(i,p)=Nph*abs(ILoad2(p).*ILoad2(p)*Rph);                                           % Watts
      % Pr= Nph*Iph^2*Rph;
        % Core Losses: W
        PCL(i,p)=(kh*fee(i)*Bmax^1.9+(ke*fee(i).^2*Bmax^2)+kc*fee(i).^1.5*Bmax^1.5)*Vst;           % Watts
        %PCL(i,p)=((kh*fe*Bmax^1.7)+(ke*fe^2*Bmax^2)+(kc*fe^1.5*Bmax^1.5))*Vst;
        %PCL(i,p)=(kh*60*Bmax^2+ke*60.^2*Bmax^2+ke*60.^1.5*Bmax^1.5)*Vst;
        % Efficiency:
        Pout(i,p)=sqrt(3)*VT2(i,p)*ILoad2(p);                              % Watts
        Pin(i,p)=Pout(i,p)+PCU(i,p)+PCL(i,p)+0.01*Pout(i,p);               % Watts
        eta(i,p)=100*abs(Pout(i,p))/abs(Pin(i,p))   ;                      % %
        % Voltage Regulation (%):
        Reg(i,p)=100*(VT2(i,1)-VT2(i,p))/VT2(i,p);                         % %
    end
end



map=[0 0 0.1
    0 0 0.15
    0 0 0.2
    0 0 0.25
    0 0 0.3
    0 0 0.35
    0 0 0.4    
    0 0 0.45
    0 0 0.5
    0 0 0.55
    0 0 0.6
    0 0 0.65
    0 0 0.7
    0 0 0.75
    0 0 0.8
    0 0 0.85
    0 0 0.9];

map2=[0 0 1
    0 0 0.9
    0 0 0.7
    0 0 0.5
    0 0 0.4
    0 0 0.3
    0 0 0.2
    0 0 0.1
    0 0 0.05];

map3=[
    0 0 0.6
    0 0 0.55
    0 0 0.5
    0 0 0.45
    0 0 0.4
    0 0 0.35
    0 0 0.3
    0 0 0.25
    0 0 0.2
    0 0 0.15
    0 0 0.1
    0 0 0.07
    0 0 0.05
    0 0 0.01];

map4=[0.1 0.1 0.1
    0.25 0.25 0.25
    0.4 0.4 0.4
    0.55 0.55 0.55
    0.7 0.7 0.7
    0.9 0.9 0.9];



map5=[0.1 0.1 0.1
    0.15 0.15 0.15
    0.2 0.2 0.2
    0.25 0.25 0.25
    0.30 0.30 0.30
    0.35 0.35 0.35
    0.40 0.40 0.40
    0.45 0.45 0.45
    0.50 0.50 0.50
    0.55 0.55 0.55
    0.60 0.60 0.60
    0.65 0.65 0.65
    0.7 0.7 0.7
    0.75 0.75 0.75
    0.80 0.80 0.80
    0.85 0.85 0.85
    0.9 0.9 0.9];


map6=[0 0.1 0.1
    0 0.15 0.15
    0 0.2 0.2
    0 0.25 0.25
    0 0.30 0.30
    0 0.35 0.35
    0 0.40 0.40
    0 0.45 0.45
    0 0.50 0.50
    0 0.55 0.55
    0 0.60 0.60
    0 0.65 0.65
    0 0.7 0.7
 ];


% Electrical frequency (Hz):
fe=Nm*wm/(4*pi);              % Hz
NPCL=120*fe/Nm;

% Eficiency:
% Copper losses:
% P2=abs(I.*I*Rph);
% Core Losses:
% Pcl2= (kh*60*Bmax^2+ke*60.^2*Bmax^2+ke*60.^1.5*Bmax^1.5)*Vst;
% Pcl= (kh*fe*Bmax^2+ke*fe.^2*Bmax^2+ke*fe.^1.5*Bmax^1.5)*Vst;

% PLOT CONFIGURATION:
figure(1)
plot(n,Erms1, 'LineStyle','-','LineWidth',3,'color',[0 1 0]);
hold on;
plot(n,abs(VT1),'LineStyle','-','LineWidth',3,'color',[0 0 1]);
hold off
xlim([0 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('V_T (V)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
%h=legend('No Load','Full-Load','Location','Northwest');
h=legend('Voltaje en vacío','Voltaje con carga','Location','Northwest');
%h.Font Size = 15;
lgd.FontWeight = 'normal';
%set(h,'LineStyle','-','LineWidth',2);
grid

figure(2)
Erms=ones(1,length(I))*Erms60;
% v2= [202.301 217.90 224.61 227.26 234.86 240.47 243.91 246.21 247.82 249.03 249.94];
% x2=100*[12.20 10.88 10.02 9.62 8.22 6.89 5.91 5.16 4.57 4.1 3.72]/10;
plot(LP,(Erms), 'LineStyle','-.','LineWidth',3,'color',[0 1 0]);
hold on;
plot(LP,(VT), 'LineStyle','-','LineWidth',3,'color',[0 0 1]);
hold off;
% hold on;
% plot(x2,v2, 'LineStyle','-','LineWidth',3,'color',[0 0 1]);
% hold off;
xlim([0 100]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('Voltaje total (V)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
%h=legend('No Load','Load','Location','south');
h=legend('Voltaje en vacío','Voltaje con carga','Location','south');
h.FontSize = 15;
lgd.FontWeight = 'normal';
%set(h,'LineStyle','-','LineWidth',2);
grid

% Copper Losses:
figure(3)
%colormap Hot;
colormap (map6)
%colormap(map)
%colormap Summer;
%colormap winter
%colormap bone
%colormap copper
%colormap pink
%colormap gray
%colormap spring
%colormap default;
%colormap parula;
%colormap([1,0,0; 1,1,1])
%colormap(flipud(colormap))
surf((LP), n, (VT2),'FaceLighting','gouraud','LineWidth',0.5);view(50,20);
set(gca,'fontname','tahoma','fontsize',15,'fontweight','light','LineWidth',2);

ylabel('S_{r} (rpm)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Carga (%)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
zlabel('V_{L}(V)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
alpha(.9)
%Voltaje con carga 

figure(4)
%colormap Hot;
colormap (map6)
%colormap(map)
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
%colormap(flipud(colormap))
xlim([0 500]);
surf(LP, n, (PCU),'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
set(gca,'fontname','tahoma','fontsize',15,'fontweight','light','LineWidth',2);
ylabel('S_{r} (rpm)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Carga (%)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
zlabel('P_{CU}(W)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
alpha(.9)
%Pérdidas del cobre 
figure(5)
%colormap Hot;

colormap (map6)
%colormap bone
%colormap(map)
%colormap Summer;
%colormap default;
%colormap([1,0,0; 1,1,1])
%colormap(flipud(colormap))
xlim([0 500]);
surf(LP, n, (PCL),'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
set(gca,'fontname','tahoma','fontsize',15,'fontweight','light','LineWidth',2);
ylabel('n (rpm)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Carga (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
zlabel('P_{CL}(W)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
alpha(.9)
%Pérdidas del núcleo 

% Efficiency:
figure(6)
%colormap Hot;
%colormap parula;
%colormap (map2);
colormap (map6)
%colormap colorcube
%colormap default;
%colormap([1,0,0; 1,1,1])
colormap(flipud(colormap))
xlim([0 500]);
zlim([20 100]);
axis tight
surf(LP, n, eta,'FaceLighting','gouraud','LineWidth',0.5);view(-50,20);
set(gca,'fontname','tahoma','fontsize',15,'fontweight','light','LineWidth',2);
ylabel('S_{r} (rpm)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Carga (%)','fontsize',24,'fontname','normal','rotation',5,'FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
zlabel('Eficiencia (%)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
alpha(.9)

% Mapa de Contorno de la eficiencia:
figure(7)
colormap (map6)
%colormap (map)
%colormap GRAY;
% colormap hot;
%colormap parula
%colormap colorcube
contourf(LP,n,eta,[0 60 80 90 95 100])
xlim([0 100]);
ylim([0 500]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('S_{r} (rpm)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
xlabel('Carga (%)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
%xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
%'HorizontalAlignment','center');
text(60,250,'\bf90% <\eta \leq 95%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(47,75,'\bf80% <\eta \leq 90%\rm','FontSize',15,'BackgroundColor','white');
text(10,80,'\bf60% <\eta \leq 80%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(1,40,'\bf\eta \leq 60%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
text(78,372,'\bf\eta > 95%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%contour(DP,(ETA),MGL,30);
%contour3(DP,(ETA),MGL,30);
alpha(.9)

% Mapa de Contorno de la Regulación de Voltaje:
figure(8)
%colormap GRAY;
colormap(map6)
% colormap hot;
contourf(LP,n,Reg,[0 1 2 5 10 50 100]);
xlim([0 100]);
ylim([25 400]);
set(gca,'fontname','tahoma','fontsize',18,'fontweight','light','LineWidth',2);
ylabel('S_{r} (rpm)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal','rotation',90,...
 'HorizontalAlignment','center');
% ylim([-150 200]);
%xlabel('Power (%)','fontsize',18,'fontname','normal','FontAngle','normal','fontweight','normal',...
%'HorizontalAlignment','center');
xlabel('Carga (%)','fontsize',24,'fontname','normal','FontAngle','normal','fontweight','normal',...
'HorizontalAlignment','center');
text(60,255,'\bf5%< V_r_e_g \leq10%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(47,145,'\bf2% < V_r_e_g \leq 5%\rm','FontSize',15,'BackgroundColor','white');
text(40,65,'\bf1% < V_r_e_g \leq 2%\rm','fontname','normal','FontSize',15,'BackgroundColor','white');
text(5,60,'\bfV_r_e_g \leq 1%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
text(76,370,'\bfV_r_e_g > 10%\rm','FontSize',15,'fontname','normal','BackgroundColor','white');
%contour(DP,(ETA),MGL,30);
%contour3(DP,(ETA),MG
alpha(.9)

figure(20)
am=xlsread('am3.csv');
Y1=bar(am(:,1),am(:,2),'FaceColor',[0 0.7 0 ],'EdgeColor',[0 0 0],'LineWidth',2);
%xlim([57 88]);
%ylim([0 40]);
%set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
%ylabel('THD_{V} [%]','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
 %   'HorizontalAlignment','center');
%xlabel('\alpha_m (%)','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
  %  'HorizontalAlignment','center');
%h = legend('THD_{V}','Location','northeast','Orientation','horizontal');
%set(h,'LineWidth',2,'fontsize',14);

hold on
%figure(21)
am=xlsread('am3.csv');
Y1=bar(am(:,1),am(:,3),'FaceColor',[0 0 0.7],'EdgeColor',[0 0 0],'LineWidth',2);
xlim([57 88]);
ylim([0 40]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('THD_{V} [%]','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('\alpha_m (%)','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('THD_{Vacío}','THD_{Carga}','Location','northeast','Orientation','horizontal');
set(h,'LineWidth',2,'fontsize',14);
grid;

figure (10)
w=[Wst,Wm,Wcu,Wrot]
wtt=Wst+Wm+Wcu+Wrot
pie3(w)
set(gca,'fontname','tahoma','fontsize',20,'fontweight','light','LineWidth',2);
h = legend('Masa de estatores','Masa de imanes','Masa de conductores', 'Masa de rotor','Location','northeast','Orientation','vertical');
set(h,'LineWidth',2,'fontsize',14);
grid;

figure (11)
Pc=[Pcl,Pr,Ps]
Pct=Pcl+Pr+Ps;
pie3(Pc)
set(gca,'fontname','tahoma','fontsize',20,'fontweight','light','LineWidth',2);
h = legend('Pérdidas del núcleo','Pérdidas del conductor','Pérdidas por fricción','Location','northeast','Orientation','vertical');
set(h,'LineWidth',2,'fontsize',14);
grid;

figure(12)
vl=xlsread('voltaje.csv');

plot(vl(:,1),vl(:,3) ,'LineStyle','-.','LineWidth',3,'color',[0 1 0]);
hold on
plot(vl(:,1),vl(:,4), 'LineStyle','-.','LineWidth',3,'color',[0 0 1]);
plot(vl(:,1),vl(:,5), 'LineStyle','-.','LineWidth',3,'color',[0 0.5 1]);
plot(n,Erms1, 'LineStyle','-','LineWidth',3,'color',[0.2 0.3 0]);
%plot(vl(:,1),vl(:,2), 'LineStyle','-','LineWidth',3,'color',[1 0 0]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Voltaje (V)','fontsize',24,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('S_{r}(rpm)','fontsize',24,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('VAB_{LL}RMS','VBC_{LL}RMS','VCA_{LL}RMS','V_{LL}RMS MATLAB','Location','Northwest','Orientation','vertical');
grid

figure(13)

%Y1=bar(v1,'stacked','FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);

%Y1=bar(v2,v3,'FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);
%Y2=bar(vl(:,1),vl(:,3),'stacked','FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);
Y4=bar(vl(:,1),vl(:,6),'stacked','FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);
hold on
Y1=bar(vl(:,1),vl(:,3),'stacked','FaceColor',[0 0 1],'EdgeColor',[0 0 0],'LineWidth',2);
xlim([0 550]);
ylim([0 1000]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Voltaje (V)','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('S_{r} (rpm)','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('VA_{MAX}','VA_{RMS}','Location','northeast','Orientation','horizontal');
set(h,'LineWidth',2,'fontsize',14);
grid;
amm=xlsread('amm2.csv');

figure(14)

%Y1=bar(v1,'stacked','FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);

%Y1=bar(v2,v3,'FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);
%Y2=bar(vl(:,1),vl(:,3),'stacked','FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);
Y6=bar(amm(:,1),amm(:,3),'stacked','FaceColor',[0 1 0],'EdgeColor',[0 0 0],'LineWidth',2);
hold on
Y7=bar(amm(:,1),amm(:,6),'stacked','FaceColor',[0 0 1],'EdgeColor',[0 0 0],'LineWidth',2);
Y8=bar(amm(:,1),amm(:,12),'stacked','FaceColor',[1 0 0],'EdgeColor',[0 0 0],'LineWidth',2);

xlim([52 98]);
ylim([0 500]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Voltaje (V)','fontsize',24,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('S_{r} (rpm)','fontsize',24,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('VA_{MAX}','VA_{RMS}','THD','Location','northeast','Orientation','horizontal');
set(h,'LineWidth',2,'fontsize',14);
grid;



figure(16)
v2=xlsread('voltajes.csv');
y1=v2(:,2.)-v2(:,3.);
y2=v2(:,2);
y5=v2(:,2.)-v2(:,4.);
y6=v2(:,2.)-v2(:,5.);
y3=v2(:,1);
y4=v2(:,5);
% y1=[v2(:,1),v2(:,3)];
% y2=[v2(:,1),v2(:,4)];
% y3=[v2(:,1),v2(:,5)];
plot(v2(:,1),v2(:,3) ,'LineStyle','-.','LineWidth',3,'color',[0 1 0]);
hold on
plot(v2(:,1),v2(:,4), 'LineStyle','-.','LineWidth',3,'color',[0 0 1]);
plot(v2(:,1),v2(:,5), 'LineStyle','-.','LineWidth',3,'color',[0 0.5 1]);
plot(n,abs(VT1),'LineStyle','-','LineWidth',3,'color',[1 0 0]);
plot(n,Erms1, 'LineStyle','-','LineWidth',3,'color',[0.2 0.3 0]);
errorbar(y3,y2,y1);
errorbar(y3,y2,y5);
errorbar(y3,y2,y6);
%plot(v2(:,1),v2(:,2), 'LineStyle','-','LineWidth',3,'color',[1 0 0]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Voltaje (V)','fontsize',24,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('S_{r}(rpm)','fontsize',24,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('VAB_{LL}RMS','VBC_{LL}RMS','VCA_{LL}RMS','V_{LLC}RMS MATLAB','V_{LLV}RMS MATLAB','Location','Northwest','Orientation','vertical');
grid

