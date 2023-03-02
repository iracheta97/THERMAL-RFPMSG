%{
  Programa para el diseño de devanados de un RF-PMSG
  Proyecto de Tesis: Eduardo Ortiz Garcia
  Algoritmo de Devanados
  V_Final
%}

clc; 
clear;
%--------------------------------------------------------------------------
% Initial Parameters
%--------------------------------------------------------------------------
% Number of slots
Ns=96;  
% Number of poles
Nm=32;                           
% Number of phases
Nph=3;                            
% Number of layers
Nlay=1;  

%--------------------------------------------------------------------------
% Evaluation of initial parameters
%--------------------------------------------------------------------------
% Number of phasors having equal phases
t=gcd(Ns,(Nm*0.5)); 

%--------------------------------------------------------------------------
% Balanced winding or Unbalanced winding
%--------------------------------------------------------------------------

%{
An unbalanced winding has a combination of number of poles and number 
of slots that does not allow to arrange the coils in such a way that 
they produce a symmetrical system of equally time-phase displaced emf's 
of identical magnitude, frequency, and waveform.
%}

% Number of winding layers
if Nlay==1                                % Single-layer winding
    y1=Ns/(2*Nph);                   
    sc=(2*y1)/t;
    ac=sc-floor(sc);
    if ac==0
         l='';
         l1=0;
         else 
         l='No se puede calcular';
         l1=1;
    end
    disp(l);
else                                      % Double-layer winding
    y1=Ns/Nph;                    
    dc=(y1/t);
    dd=dc-floor(dc);
    if dd==0
        ll='';
        l1=0;
        else 
        ll='No se puede calcular';
        l1=1;
    end   
    disp(ll);  
end

%--------------------------------------------------------------------------
if l1==0  % Symmetry conditions
%--------------------------------------------------------------------------
% Type of winding
wtype=1;
% Number of empty slots 
Nes=0;                             
% Number of slots per pole - coil span 
Nsm=Ns/Nm;                
% Diamond-coil winding is chosen, with shorted pitch
Yq=floor(Nsm);
% Number of slots per pole per phase:
Nspp=Ns/(Nph*Nm);
% Number of slots per phase:
Nsp=Ns/Nph;
% Number of coiled windings:
Nwc=(wtype*(Ns-Nes)/(2*Nph));
e=Nwc-floor(Nwc);
if e~=0
    Nes=Ns-2*Nph*floor(Nwc)/wtype;          % Number of empty slots
end
% Number of wound slots per pole per phase:
q=(Ns-Nes)/(Nph*Nm);
a=floor(q);
z=wtype*0.5*Nm*(q-a);
% Number of repetitions:
t1=gcd(round(z),0.5*Nm);
pp= Nm/(2*t1);
repeticiones=t1;
if t1==1                                                   
    g= wtype*Ns/(2*Nph*t);        % Normal - Non- Reduced Systems 
else                              
    g= Ns/(2*Nph*t);              % Reduced subsystems
end
f=g-floor(g);
if f==0;                                              
    sym='Devanado Simetrico';       % Symmetrical Winding 
else
    sym='Devanado Asimetrico';      % Assymetrical Winding
end

if Nspp == q                      % Selecciona el tipo de devanado
    A='Devanado de Ranura Integral';
elseif Nspp <= 1
    A='Devanado Concentrado de Ranura Fraccional';
elseif Nspp >= 1
    A='Devanado Distribuido de Ranura Fraccional';
end
N1=Ns/t1;  
%--------------------------------------------------------------------------
% Factor de devando 
%--------------------------------------------------------------------------
al=(pi/(Nm*Nspp));
alfa2=((0.5*Nm)*(2*pi))/Ns;
ed=Nsm-floor(Nsm);
% Factor de devanado 
Kw=[]; Kw1=[]; Kw2=[];
for r=1:2:20
    Kd=(sin((r*Nspp)*(alfa2/2)))/(Nspp*sin(r*alfa2/2));
    %Kd=sin(pi/(2*Nph))/(Nspp*sin(pi/(2*Nph*Nspp)));
    if ed==0;
        if Nlay==1    
    Kp=sin((r*(Nsm)*pi)/((Nsm)*2));  % for v odd
        else
    Kp=sin((r*(dc)*pi)/((Nsm)*2));  % for v odd
        end
   
    else 
        if Nlay==1    
    Kp=sin((r*(Nsm)*pi)/((Nsm)*2));  % for v odd
        else
    Kp=sin((r*(dc)*pi)/((Nsm)*2));  % for v odd
        end
    
    end
    KW=abs(Kd*Kp);
    Kw=[Kw KW];
    KW11=Kw';
    KW1=abs(Kd*Kp)/r;   
    Kw1=[Kw1 KW1];
    KW12=Kw1';
    KW2=(Kd*Kp)/r;   
    Kw2=[Kw2 KW2];
    KW13=Kw2'
    
end
figure(1)
Y1=bar(KW11,'FaceColor',[0 .5 .5],'EdgeColor',[0 0 0],'LineWidth',1,'BarWidth',0.8);
% bar(Kw)
title('Factor de devanado', 'fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
ylabel('Amplitud','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('Orden de armónico','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlim([1 10])
grid on
figure(2)
Y2=bar(KW12,'FaceColor',[0 .5 .5],'EdgeColor',[0 0 0],'LineWidth',1,'BarWidth',0.8);
%bar(Kw1)
title('Fase única - MMF', 'fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
ylabel('Amplitud','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('Orden de armónico','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
xlim([1 10])
grid on
figure(3)
Y3=bar(KW13,'FaceColor',[0 .5 .5],'EdgeColor',[0 0 0],'LineWidth',1,'BarWidth',0.8);
% bar(Kw2)
title('Tres fases - MMF', 'fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
ylabel('Amplitud','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal','rotation',90,...
    'HorizontalAlignment','center');
xlabel('Orden de armónico','fontsize',18,'fontname','times','FontAngle','normal','fontweight','normal',...
    'HorizontalAlignment','center');
ylim([0 1])
xlim([1 10])
grid on
%--------------------------------------------------------------------------
% Tabla de distribución de devanados (TDD) 
%--------------------------------------------------------------------------
Fasor=Ns/t;                       % Número de de fasores compuestos
alfa=t*360/Ns;                    % Ángulo electrico alfa 
alfa2=((0.5*Nm)*(2*pi))/Ns;       % Ángulo de la ranura

m2=rem(Nph,2);
if m2==0
    shup=(Nph)/2-1;        % Sift up
else
    shup=(Nph-1)/2;        % Sift up
end

R=Fasor*(1:t)-Fasor;
for n=1:Fasor
    Table(1,n)=alfa*(n-1);
    Table((1:t)+1,n)=R+n;
end

M=zeros(t1,Ns/t1);

cte=0;
for q=1:t1
for i=1:Ns/t1
    ct(i)=cte;
    c1=Ns/t1;
    ind=mod((q-1)*c1+(0.5*Nm/t1)*(i-1)+1,q*Ns/t1);
  
    if ind==0
         ind=q*Ns/t1;
    end
    % Pointer:
    q1=ind-(q-1)*Ns/t1;
    %x(i)=M(q,ind)
    if M(q,q1)~=0      
            ind=ind+cte;
    end
    if rem(i,Ns/t)==0
        cte=cte+1;
    end
    % Update Pointer:
    q1=ind-(q-1)*Ns/t1;
    M(q,q1)=i+(q-1)*Ns/t1;
    y(:,i)=[i;ind];
    i=i+1;
    ind=ind+1;   
end
end

M1=zeros(Nph,2/wtype*Nwc/t1);
Nc1=floor(Ns/(Nph*t1)/2);
Nc2=round(Ns/(Nph*t1)/2);
Pos=1:Nc2;
Neg=Nc2+1:Ns/(Nph*t1);
WDT=zeros(Nph,Ns/(Nph));
Nsp2=Ns/(2*Nph);


Ma=[];

for q=1:t1
for i=1:Nph
     M1(i,:,q)=M(q,2/wtype*Nwc/t1*(i-1)+1:2/wtype*Nwc/t1*i)
     M1(i,Neg,q)=-M1(i,Neg,q)
end
     vec=Pos+length(Pos)*(q-1);
     WDT(:,vec)=M1(:,Pos,q);
     x1=length(Neg)*(q-1)+1:length(Neg);
     vec2= t1*length(Pos)+[1:length(Neg)]+length(Neg)*(q-1);
     WDT(:,vec2)=M1(:,Neg,q);
end

shift=[shup+1:Nph 1:shup];
w=WDT;
WDT(:,length(Pos)*(t1)+1:Nsp)=...
    WDT(shift,length(Pos)*(t1)+1:Nsp);

if Nspp==1

[r1 c1]=find(WDT<0);
R1=length(r1); C1=length(c1);
WDT3=sort(abs(WDT),2);
WDT4=WDT3;
for i=1:R1
    dd=abs(WDT(r1(i),c1(i)));
    [r11 c11]=find(WDT3==dd);
    WDT3(r11,c11)=-WDT3(r11,c11);
end

W1=[];W2=[];A=0;B=0;C=0;
for i=1:Ns
    % Find locations of M1:
    [s t1]=find(WDT3==i|WDT3==-i);
    %
    if s==1                   % Phase A
           A=A+1;
           if WDT3(s,t1)<0
           W1=strcat(W1,'X');    
           else
           W1=strcat(W1,'A');
           end
           elseif s==2        % Phase B
           B=B+1;
           if WDT3(s,t1)<0
           W1=strcat(W1,'Y');
           else
           W1=strcat(W1,'B');
           end
           elseif s==3        % Phase C
           C=C+1;
           if WDT3(s,t1)<0
           W1=strcat(W1,'Z');
           else
           W1=strcat(W1,'C');
           end
    end
end
    
WW=strcat(W1,W2);
 
[r1 c1]=find(WDT<0);
R1=length(r1);
C1=length(c1);
WDT2=sort(abs(WDT),2);
WDT3=WDT2;
for i=1:R1
    dd=abs(WDT(r1(i),c1(i)));
    [r11 c11]=find(WDT2==dd);
    WDT2(r11,c11)=-WDT2(r11,c11);

end

WDT4=WDT3+3;
n=rem((WDT4),Ns);
n1=n;
n1(n1==0)=Ns;

DT=-ones(Nph,Ns/(Nph));
D=WDT2.*DT;
D1=abs(D);
T=(n1.*D);
WDT5=T./D1;
  
WDT6=[WDT2,WDT5];
D=[WDT6(:)]';
d=(Ns/(Nph));

WDT7=[];
for j=1:d
s=[WDT2(:,j) WDT5(:,j)];
WDT7=[WDT7 s];      
end
     
W1=[];W2=[];A=0;B=0;C=0;
for i=1:Ns
    % Find locations of M1:
    [s t1]=find(WDT2==i|WDT2==-i);
    %
    if s==1                   % Phase A
           A=A+1;
           if WDT2(s,t1)<0
           W1=strcat(W1,'X'); 
           else
           W1=strcat(W1,'A');
           end       
           elseif s==2        % Phase B
           B=B+1;
           if WDT2(s,t1)<0
           W1=strcat(W1,'Y');
           else
           W1=strcat(W1,'B');
           end
           elseif s==3        % Phase C
           C=C+1;
           if WDT2(s,t1)<0
           W1=strcat(W1,'Z');
           else
           W1=strcat(W1,'C');
           end
    end
end
    
WW=strcat(W1,W2);

W3=[];W4=[];A=0;B=0;C=0;
for i=1:Ns
    % Find locations of M1:
    [s t1]=find(WDT5==i|WDT5==-i);
    %
    if s==1                   % Phase A
           A=A+1;
           if WDT5(s,t1)<0
           W3=strcat(W3,'X'); 
           else
           W3=strcat(W3,'A');
           end
           elseif s==2        % Phase B
           B=B+1;
           if WDT5(s,t1)<0
           W3=strcat(W3,'Y');
           else
           W3=strcat(W3,'B');
           end
           elseif s==3        % Phase C
           C=C+1;
           if WDT5(s,t1)<0
           W3=strcat(W3,'Z');
           else
           W3=strcat(W3,'C');
           end
    end
end
    
WW1=strcat(W3,W4);

WW3=[];
for k=1:Ns
R=[WW(:,k) WW1(:,k)];
WW3=[WW3 R];    
end

else
    
z=(length(WDT)/2);
zx=length(WDT);

WDT7=[];
for j=1:z
s=[WDT(:,j)];
WDT7=[WDT7 s]; 
end

WDT8=[];
for p=z+1:zx
s=[WDT(:,p)];
WDT8=[WDT8 s];  
end

WDT711=[];
for j=2:2:z
s=[WDT7(:,j)];
WDT711=[WDT711 s]; 
end

WDT712=[];
for j=1:2:z
s=[WDT7(:,j)];
WDT712=[WDT712 s]; 
end

B1=WDT711(1,:);
WDT711(1,:)=WDT711(2,:);
WDT711(4,:)=B1;
WDT711(1,:)=[] ;
WDT711=WDT711*-1;

WDT811=[];
for j=3:2:z
s=[WDT8(:,j)];
WDT811=[WDT811 s]; 
end

WDT812=[];
for j=2:2:z
s=[WDT8(:,j)];
WDT812=[WDT812 s]; 
end

B2=WDT811(1,:);
WDT811(1,:)=WDT811(2,:);
WDT811(4,:)=B2;
WDT811(1,:)=[] ;
WDT811=WDT811*-1;

WDTT=[WDT712,WDT711];
WDTT1=[WDT812,WDT811];
WDTT2=[WDTT,WDTT1];
FIN=[WDTT2 WDT8(:,1)];

[r1 c1]=find(FIN<0);
R1=length(r1); C1=length(c1);
WDT3=sort(abs(FIN),2);
WDT4=WDT3;
for i=1:R1
    dd=abs(FIN(r1(i),c1(i)));
    [r11 c11]=find(WDT3==dd);
    WDT3(r11,c11)=-WDT3(r11,c11);
end

W1=[];W2=[];A=0;B=0;C=0;
for i=1:Ns
    % Find locations of M1:
    [s t1]=find(WDT3==i|WDT3==-i);
    %
    if s==1                   % Phase A
           A=A+1;
           if WDT3(s,t1)<0
           W1=strcat(W1,'X');    
           else
           W1=strcat(W1,'A');
           end
           elseif s==2        % Phase B
           B=B+1;
           if WDT3(s,t1)<0
           W1=strcat(W1,'Y');
           else
           W1=strcat(W1,'B');
           end
           elseif s==3        % Phase C
           C=C+1;
           if WDT3(s,t1)<0
           W1=strcat(W1,'Z');
           else
           W1=strcat(W1,'C');
           end
    end
end
    
WW=strcat(W1,W2);

[r1 c1]=find(FIN<0);
R1=length(r1);
C1=length(c1);
WDT2=sort(abs(FIN),2);
WDT3=WDT2;
for i=1:R1
    dd=abs(FIN(r1(i),c1(i)));
    [r11 c11]=find(WDT2==dd);
    WDT2(r11,c11)=-WDT2(r11,c11);

end

WDT4=WDT3+3;
n=rem((WDT4),Ns);
n1=n;
n1(n1==0)=Ns;

DT=-ones(Nph,Ns/(Nph));
D=WDT2.*DT;
D1=abs(D);
T=(n1.*D);
WDT5=T./D1;
  
WDT6=[WDT2,WDT5];
D=[WDT6(:)]';
d=(Ns/(Nph));

WDT7=[];
for j=1:d
s=[WDT2(:,j) WDT5(:,j)];
WDT7=[WDT7 s];      
end
     
W1=[];W2=[];A=0;B=0;C=0;
for i=1:Ns
    % Find locations of M1:
    [s t1]=find(WDT2==i|WDT2==-i);
    %
    if s==1                   % Phase A
           A=A+1;
           if WDT2(s,t1)<0
           W1=strcat(W1,'X'); 
           else
           W1=strcat(W1,'A');
           end       
           elseif s==2        % Phase B
           B=B+1;
           if WDT2(s,t1)<0
           W1=strcat(W1,'Y');
           else
           W1=strcat(W1,'B');
           end
           elseif s==3        % Phase C
           C=C+1;
           if WDT2(s,t1)<0
           W1=strcat(W1,'Z');
           else
           W1=strcat(W1,'C');
           end
    end
end
    
WW=strcat(W1,W2);

W3=[];W4=[];A=0;B=0;C=0;
for i=1:Ns
    % Find locations of M1:
    [s t1]=find(WDT5==i|WDT5==-i);
    %
    if s==1                   % Phase A
           A=A+1;
           if WDT5(s,t1)<0
           W3=strcat(W3,'X'); 
           else
           W3=strcat(W3,'A');
           end
           elseif s==2        % Phase B
           B=B+1;
           if WDT5(s,t1)<0
           W3=strcat(W3,'Y');
           else
           W3=strcat(W3,'B');
           end
           elseif s==3        % Phase C
           C=C+1;
           if WDT5(s,t1)<0
           W3=strcat(W3,'Z');
           else
           W3=strcat(W3,'C');
           end
    end
end
   
WW1=strcat(W3,W4);

WW3=[];
for k=1:Ns
R=[WW(:,k) WW1(:,k)];
WW3=[WW3 R];    
end

end

else
    disp('Cambia los parametros de entrada');
    end
















