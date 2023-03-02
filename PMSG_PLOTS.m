function PMSG_PLOTS(PMSG)
% Help: This function plots thte performance features of the RF_PMSG
% Created by: Dr. Reynaldo Iracheta Cortez
% Date: December 10th, 2020
% Inputs:
% PMSG: Structure Data of the PMSG
% Outputs:
% PLOT of temperature as function of load percentage
%--------------------------------------------------------------------------
% Load percentage:
%--------------------------------------------------------------------------
PL=0:1e-3:1;                                         %
Iph=PMSG.ELECTRICAL.Iph;                             % A
TF=[];
TS=[];
% Initialization of vector of currents:
for i=1:length(PL)  
        PMSG.ELECTRICAL.Iphx=Iph*PL(i);              % A
        [PMSG0]=PMSG_THERMAL(PMSG);
        % Tempertatures obtained with the full model:
        TF(:,i)=PMSG0.THERMAL.TF;
        % Temperatures obtained with the simplified model:
        TS(:,i)=PMSG0.THERMAL.TS;                     % °C
end
%--------------------------------------------------------------------------
% Environmental temperature: (°C) Full load
%--------------------------------------------------------------------------
TT=0.1:1:50;                                          % °C
TF1=[];
TS1=[];
for i=1:length(TT)
    PMSG.ELECTRICAL.Iph=Iph;                          % A
    PMSG.THERMAL.Tair=TT(i);                          % °C
    [PMSG1]=PMSG_THERMAL(PMSG);
    % Temperatures obtained with the full model:
    TF1(:,i)=PMSG1.THERMAL.TF;
    % Temperatures obtained with the simplified model:
    TS1(:,i)=PMSG1.THERMAL.TS;                     % °C
end
%--------------------------------------------------------------------------
% Temperatures Full Model:
%--------------------------------------------------------------------------
% Entrehierro: Air Gap
figure(1)
PL1=100*PL;
plot(PL1,TF(13,:),'color',[0 0 0],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
hold on;
% Yugo del estator:
%plot(PL1,TF(2,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(3,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
%plot(PL1,TF(4,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(5,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Diente:
plot(PL1,TF(6,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(7,:),'color',[0 0 0],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% plot(PL1,TF(8,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(PL1,TF(9,:),'color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% Zapato:
%plot(PL1,TF(10,:),'color',[0 1 0],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(11,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(12,:),'color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(13,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Magnet Temperature:
plot(PL1,TF(14,:),'color',[1 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
%plot(PL1,TF(20,:),'color',[1 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Rotor yoke temperature is equal to magnet temperature:
%plot(PL1,TF(15,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(16,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% End winding:
%plot(PL1,TF(17,:),'color',[1 0 1],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(18,:),'color',[1 0 1],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Devanado Centro:
plot(PL1,TF(19,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
%plot(PL1,TF(8,:),'color',[0 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(9,:),'color',[0 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);

%plot(s,Tmech,'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
hold off;
%ylim([0 5000]);
xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Temperatura (°C)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('Carga [%]','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
% h = legend('Entrehierro, Centro','Yugo del estátor 1', 'Yugo del estátor 2','Yugo del estátor 3', 'Yugo del estátor 4',...
%            'Diente 1', 'Diente 2','Zapato','Entrehierro, Superior','Entrehierro,Inferior',...
%            'Imán 1','Imán 2','Rotor 1', 'Rotor 2','Fin de devanado 1','Fin del devanado 2',...
%            'Devanado Centro (Axial)', 'Devanado Centro (Radial 1)','Devanado Centro (Radial 2)','Location','northwest','Orientation','vertical','NumColumns',2);
h = legend('Entrehierro',...
           'Estátor',...
           'Imán-Rotor','Devanado','Location','northwest','Orientation','vertical','NumColumns',2);
%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
title(h,['Distribución de Temperaturas en un RF-PMSG de 10 kW']);
grid;
%--------------------------------------------------------------------------
% Temperatures Simplified Model:
%--------------------------------------------------------------------------
figure(2)
PL1=100*PL;
% Entrehierro:
plot(PL1,TS(8,:),'color',[0 0 0],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
hold on;
% Yugo del estátor:
%plot(PL1,TS(2,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TS(3,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% Diente:
plot(PL1,TS(4,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Entrehierro:
%plot(PL1,TS(8,:),'color',[1 0 0],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Imán:
plot(PL1,TS(9,:),'color',[1 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Devanado:
plot(PL1,TS(5,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TS(6,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Devanado 2:
% plot(PL1,TS(9,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TS(10,:),'color',[0.49 .18 .56],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Rotor
%plot(PL1,TS(10,:),'color',[1 0 0],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Devanado Axial:
%plot(PL1,TS(12,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);

%plot(s,Tmech,'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
hold off;
%ylim([0 5000]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Temperatura (°C)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('Carga [%]','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('Entrehierro',...
           'Estator','Imán-Rotor','Devanado','Location','northwest','Orientation','vertical','NumColumns',2);%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
%title(h,['Distribución de Temperaturas en la Máquina de Imanes Permanentes']);
title(h,['Distribución de Temperaturas en un RF-PMSG de 10 kW']);
grid;
%--------------------------------------------------------------------------
% Temperatures Full Model:
%--------------------------------------------------------------------------
% Entrehierro: Air Gap
figure(3)
PL1=100*PL;
plot(TT,TF1(13,:),'color',[0 0 0],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
hold on;
% Yugo del estator:
%plot(PL1,TF(2,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(3,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
%plot(PL1,TF(4,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(5,:),'color',[1 0 1],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Diente:
plot(TT,TF1(6,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(7,:),'color',[0 0 0],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% plot(PL1,TF(8,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% plot(PL1,TF(9,:),'color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% Zapato:
%plot(PL1,TF(10,:),'color',[0 1 0],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(11,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(12,:),'color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(13,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Magnet Temperature:
plot(TT,TF1(14,:),'color',[1 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
%plot(PL1,TF(20,:),'color',[1 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Rotor yoke temperature is equal to magnet temperature:
%plot(PL1,TF(15,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(16,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% End winding:
%plot(PL1,TF(17,:),'color',[1 0 1],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(18,:),'color',[1 0 1],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Devanado Centro:
plot(TT,TF1(19,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
%plot(PL1,TF(8,:),'color',[0 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TF(9,:),'color',[0 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);

%plot(s,Tmech,'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
hold off;
%ylim([0 5000]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Temperatura (°C)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('Temperatura Ambiente [°C]','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
% h = legend('Entrehierro, Centro','Yugo del estátor 1', 'Yugo del estátor 2','Yugo del estátor 3', 'Yugo del estátor 4',...
%            'Diente 1', 'Diente 2','Zapato','Entrehierro, Superior','Entrehierro,Inferior',...
%            'Imán 1','Imán 2','Rotor 1', 'Rotor 2','Fin de devanado 1','Fin del devanado 2',...
%            'Devanado Centro (Axial)', 'Devanado Centro (Radial 1)','Devanado Centro (Radial 2)','Location','northwest','Orientation','vertical','NumColumns',2);
h = legend('Entrehierro',...
           'Estátor',...
           'Imán-Rotor','Devanado','Location','northwest','Orientation','vertical','NumColumns',2);
%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
title(h,['Distribución de Temperaturas en un RF-PMSG de 10 kW']);
grid;
%--------------------------------------------------------------------------
% Temperatures Simplified Model:
%--------------------------------------------------------------------------
figure(4)
PL1=100*PL;
% Entrehierro:
plot(TT,TS1(8,:),'color',[0 0 0],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
hold on;
% Yugo del estátor:
%plot(PL1,TS(2,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TS(3,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
% Diente:
plot(TT,TS1(4,:),'color',[0 0.5 0.5],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Entrehierro:
%plot(PL1,TS(8,:),'color',[1 0 0],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Imán:
plot(TT,TS1(9,:),'color',[1 0 0],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Devanado:
plot(TT,TS1(5,:),'color',[0.49 .18 .56],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TS(6,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Devanado 2:
% plot(PL1,TS(9,:),'color',[0.49 .18 .56],'LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
%plot(PL1,TS(10,:),'color',[0.49 .18 .56],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);
% Rotor
%plot(PL1,TS(10,:),'color',[1 0 0],'LineStyle','--','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',1);
% Devanado Axial:
%plot(PL1,TS(12,:),'color',[0 .5 .5],'LineStyle','-.','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5],'MarkerSize',0.001);

%plot(s,Tmech,'color',[0.49 .18 .56],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.49 .18 .56],'MarkerSize',8);
% %plot(s,[VL23],'color',[0.8500 0.3250 0.0980],'Marker','s','LineStyle','-','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
hold off;
%ylim([0 5000]);
% %xlim([100 400]);
set(gca,'fontname','tahoma','fontsize',16,'fontweight','light','LineWidth',2);
ylabel('Temperatura (°C)','fontsize',18,'fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlabel('Temperatura ambiente [°C]','fontsize',18,'fontname','tahoma','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
h = legend('Entrehierro',...
           'Estator','Imán-Rotor','Devanado','Location','northwest','Orientation','vertical','NumColumns',2);%title(h,['Devanado-Prueba,' newline 'Análisis FEA'])
%title(h,['Distribución de Temperaturas en la Máquina de Imanes Permanentes']);
title(h,['Distribución de Temperaturas en un RF-PMSG de 10 kW']);
grid;