function[GEN]= PMSG_ELECTRICAL(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
wm=GEN.wr.Data;
Emax=GEN.Emax.Data;             % Volts
fe=GEN.fe.Data;                 % Hz
we= 2*pi*fe;                    % rad/s
%--------------------------------------------------------------------------
% Topological Constraints:
%--------------------------------------------------------------------------
Nm=GEN.CONSTRAINTS.Nm;
Nph=GEN.CONSTRAINTS.Nph;
Nspp=GEN.CONSTRAINTS.Nspp;
Ns=GEN.CONSTRAINTS.Ns;
Nlay=GEN.CONSTRAINTS.Nlay;
Nsp=Ns/Nph;
%--------------------------------------------------------------------------
% Geometry:
%--------------------------------------------------------------------------
g=GEN.GEOMETRY.g;
Rso=GEN.GEOMETRY.Rso;
Rsi=GEN.GEOMETRY.Rsi;
Rro=GEN.GEOMETRY.Rro;
Rri=GEN.GEOMETRY.Rri;
Rsb=GEN.GEOMETRY.Rsb;
lm=GEN.GEOMETRY.lm;
Tp=GEN.GEOMETRY.Tp;
Ts=GEN.GEOMETRY.Ts;
Tc=GEN.GEOMETRY.Tc;
L=GEN.GEOMETRY.L;
am=GEN.GEOMETRY.am;
Cphi=GEN.GEOMETRY.Cphi;
Tetas=GEN.GEOMETRY.Tetas;
As=PMSG.GEOMETRY.As;
d3=PMSG.GEOMETRY.d3;
d2=PMSG.GEOMETRY.d2;
d1=PMSG.GEOMETRY.d1;
ds=PMSG.GEOMETRY.ds;
wsi=PMSG.GEOMETRY.wsi;
wsb=PMSG.GEOMETRY.wsb;
wt=PMSG.GEOMETRY.wt;
ws=PMSG.SHOE.ws;
%--------------------------------------------------------------------------
% Magnet:
%--------------------------------------------------------------------------
miur=GEN.MAGNET.miuR;
kc=GEN.MAGNET.kc;
Wm=GEN.MAGNET.Wm;
miuR=GEN.MAGNET.miuR;
%--------------------------------------------------------------------------
% Steel:
%--------------------------------------------------------------------------
Kst=GEN.STEEL.Kst;
rhobi=GEN.STEEL.rhobi;          % kg/m3
%--------------------------------------------------------------------------
% Winding:
%--------------------------------------------------------------------------
kcp= GEN.WINDING.kcp;
rhocu= GEN.WINDING.rhocu;
rhocux=GEN.WINDING.rhocux;      % Ohms-m
Dcu=GEN.WINDING.Dcu;            % kg/m3
dcu= GEN.WINDING.dcu;           % mm
%--------------------------------------------------------------------------
% Shoe:
%--------------------------------------------------------------------------
ws= GEN.SHOE.ws;
%--------------------------------------------------------------------------
% Flux:
%--------------------------------------------------------------------------
Bg= GEN.FLUX.Bg;
Bgx= GEN.FLUX.Bgx;
%--------------------------------------------------------------------------
% Factors:
%--------------------------------------------------------------------------
kd= GEN.FACTORS.kd;
kp= GEN.FACTORS.kp;
ks= GEN.FACTORS.ks;
%--------------------------------------------------------------------------
% Electrical:
%--------------------------------------------------------------------------
Jmax=GEN.ELECTRICAL.Jmax;
%--------------------------------------------------------------------------
% Winding:
%--------------------------------------------------------------------------
Beta=PMSG.WINDING.Beta;
%--------------------------------------------------------------------------
% Thermal:
%--------------------------------------------------------------------------
Tx=PMSG.THERMAL.Tcu;
Txm=PMSG.THERMAL.Tm;
a=PMSG.MAGNET.a;
% -------------------------------------------------------------------------
% ELECTRICAL PARAMETERS:
%--------------------------------------------------------------------------
% Number of turns per slot:
ctex=Nm*kd*kp*ks*Bgx*L*Rro*Nspp;
GEN.ELECTRICAL.ctex=ctex;
cte=Nm*kd*kp*ks*Bg*L*Rro*Nspp;
GEN.ELECTRICAL.cte=cte;
if Nlay==1
    % Considering the temperature dependence:
    nsx=floor(Emax/(ctex*wm));
    nsx=nsx+1;
    % Without thermal dependence:
    ns=floor(Emax/(cte*wm));
    ns=ns+1;   
else 
    % Considering the temperature dependence:
    nsx=floor(Emax/(ctex*wm));
    nsx=nsx+1;
    nsx=nsx/2;
    % Without thermal dependence:
    ns=floor(Emax/(cte*wm));
    ns=ns+1;
    ns=ns/2;   
end 

if Nlay==1
GEN.ELECTRICAL.nsx=nsx;
GEN.ELECTRICAL.ns=ns;
% Emax:
Emaxx=ctex*wm*ns;
Emax=cte*wm*ns;
elseif Nlay==2
GEN.ELECTRICAL.nsx=nsx;
GEN.ELECTRICAL.ns=ns;
nsx=nsx*2;
ns=ns*2;
% Emax:
Emaxx=ctex*wm*nsx;
Emax=cte*wm*ns;
end
GEN.ELECTRICAL.Emaxx=Emaxx;
GEN.ELECTRICAL.Emax=Emax;
% Peak Slot Current:
Tin=GEN.Tem.Data;
Isx= Tin/(ctex);
Is= Tin/(cte);
GEN.ELECTRICAL.Isx=Isx;
GEN.ELECTRICAL.Is=Is;
% Conductor current (A):
Icx=Isx/nsx;
Ic=Is/ns;
GEN.ELECTRICAL.Icx=Icx;
GEN.ELECTRICAL.Ic=Ic;
% Conductor Area (mm2):
%Ac=Is/(Jmax*1e6);
Ac=Isx/Jmax/1e6;                     % m2
%GEN.ELECTRICAL.Ac=Ac/1e6;           % m2
GEN.ELECTRICAL.Ac=Ac;                % m2
if Ac<As
    kcp=Ac/As;
    GEN.FACTORS.kcp=kcp;
else
    error('Aumentar los diámetros Dso, Dsi, Dri, Dro del PMSG');
end
% One conductor Area:
Ac1=pi*dcu^2/4;                     % mm2
GEN.ELECTRICAL.Ac1=Ac1;             % Ac1
% Current allowed per single conductor:
Ic1= Jmax*Ac1;                      % Amp
GEN.ELECTRICAL.Ic1=Ic1;             % Amp
% Number of parallel conductors:
npx= Ac*1e6/(nsx*Ac1);              %
np= Ac*1e6/(ns*Ac1);                %
GEN.ELECTRICAL.npx=npx;             % Amp
GEN.ELECTRICAL.np=np;               % Amp
% Phase current:
Iphx=Isx/(Nph*nsx);                 % A
GEN.ELECTRICAL.Iphx=Iphx;
Iph=Is/(Nph*ns);                    % A
GEN.ELECTRICAL.Iph=Iph;
% Peak conductor current density:
Jcx=Isx/(kcp*As);                   % A/m2
GEN.ELECTRICAL.Jcx=Jcx;
Jc=Is/(kcp*As);                     % A/m2
GEN.ELECTRICAL.Jc=Jc;
% Peak slot flux density:
miu0=4*pi*1e-7;     
Bsmaxx= miu0*Isx/ws;                       % H/m
GEN.ELECTRICAL.Bsmaxx=Bsmaxx;              % Wb/m2
Bsmax= miu0*Is/ws;                         % H/m
GEN.ELECTRICAL.Bsmax=Bsmax;                % Wb/m2
% Slot Resistance (Ohms):
% Conductor area (m2):
%Ac=kcp*As;
% Conductor heigth (m):
wc=(wsi+wsb)/2;
h= kcp*d3;                                  % m
% Conductor diameter (m):
%dcu=sqrt(4*Ac/(ns*pi));                    % m
% Updating rhocu:
rhocux=rhocu*(1+Beta*(Tx-20));
GEN.WINDING.rhocux=rhocux;                  % Ohms/m
% Skin depth:
gamx= sqrt(2*rhocux/(we*miu0));             % Copper profundity
gam= sqrt(2*rhocu/(we*miu0));               % Copper profundity
% DC Resistance:
Acsx=Ac/npx;                                % m2
Acs=Ac/np;                                  % m2
Rsx=rhocux*(nsx*nsx)*L/(Acs*npx);           % Ohms
Rs=rhocu*(ns*ns)*L/(Acs*np);                % Ohms
% AC Resistance:
% Estas fórmulas no funcionan bien:
%Rac1= (ds*rhocu*L*dcu^2*ns^2/(9*gam^4*wsb))
%deltae=(1/9)*(h^2/gam^2)*(dcu^2/gam^2)      %
Racx= rhocux*(nsx)*L/(pi*dcu*gamx*npx);      % Ohms
Rac= rhocu*(ns)*L/(pi*dcu*gam*np);           % Ohms
% Total Resistance:
Rssx=sqrt(Rsx^2+Racx^2);                      % Ohms
Rss=sqrt(Rs^2+Rac^2);                         % Ohms
%Rss1=sqrt(Rs^2+Rac1^2)                       % Ohms
GEN.ELECTRICAL.Rs=Rss;                         % Ohms
GEN.ELECTRICAL.Rsx=Rssx;                       % Ohms
% End-turn Resistance (Ohms):
% DC:
Rex=rhocux*(nsx*nsx)*pi*Tc/(2*Acsx*npx);      % Ohms
Re=rhocu*(ns*ns)*pi*Tc/(2*Acs*np);            % Ohms
% AC:
Re_acx= rhocux*nsx*pi*Tc/(pi*dcu*gamx*npx);   % Ohms
Re_ac= rhocu*ns*pi*Tc/(pi*dcu*gam*np);        % Ohms
% Total resistance:
Reex=sqrt(Rex^2+Re_acx^2);                    % Ohms
Ree=sqrt(Re^2+Re_ac^2);                       % Ohms
GEN.ELECTRICAL.Re=Ree;                        % Ohms
GEN.ELECTRICAL.Rex=Reex;                      % Ohms
% Phase Resistance (Ohms):
Rphx=Nsp*(Rssx+Reex);                        % Ohms
Rph=Nsp*(Rss+Ree);                           % Ohms
GEN.ELECTRICAL.Rphx=Rphx;                    % Ohms
GEN.ELECTRICAL.Rph=Rph;                      % Ohms
% Phase Resistance (75°C):
% Resistance at 75°C:
rhocu75=GEN.WINDING.rhocu75;                % Ohms-m
Rph75=Rph*rhocu75/rhocu;
GEN.ELECTRICAL.Rph75=Rph75;                 % Ohms
%--------------------------------------------------------------------------
% Inductance:
%--------------------------------------------------------------------------
% Air-Gap Inductance:
kc=GEN.MAGNET.kc;
gc=g+lm/miuR;
GEN.MAGNET.gc=gc;
gef= kc*gc;
GEN.MAGNET.gef=gef;
% Air-Gap Inductance (Inductance):
Lg= (ns*ns)*miur*miu0*L*Tc*kd/(4*(lm+miur*kc*g));      % H
Lgx= (nsx*nsx)*miur*miu0*L*Tc*kd/(4*(lm+miur*kc*g));   % H
GEN.ELECTRICAL.Lg=Nsp*Lg; 
GEN.ELECTRICAL.Lgx=Nsp*Lgx;  
% Slot Leakage Inductance:
% Ls= (ns*ns)*[miu0*d3^2*L/(3*(As))+ 2*miu0*d2*L/(ws+wsi) + miu0*d1*L/ws];  % H
% Xls= Nsp*we*Ls;
% Paths 5, 6 and 7:
dw=0.765*d3;
P1=miu0*gef*L/(wt+ws);
P2=miu0*d1*L/(ws);
P3=2*miu0*d2*L/(ws+wsi);
P4=miu0*(d3-dw)*L/(wsi);
P44=miu0*dw*L/(3*(wsi+wsb));
P5=miu0*(d3/2)*log(1+pi*wt/wsi)/(3*pi);
P6=miu0*log(1+pi*wt/wsi)*(d3/2)/pi;
P7=miu0*log(1+pi*wt/wsi)*(d1+d2)/pi;
Ls1=(ns*ns)*[P1+P2+P3+P4+P44+2*P5+2*P6+2*P7];          % H
Ls1x=(nsx*nsx)*[P1+P2+P3+P4+P44+2*P5+2*P6+2*P7];       % H
GEN.ELECTRICAL.Ls=Nsp*Ls1;
GEN.ELECTRICAL.Lsx=Nsp*Ls1x;  
% End Leakage Inductance:
% Segment length at the end of the turn:
Lseg=2*pi*(Rsi+ds-0.5*dw)/Ns;                       % m
% dr=d3; lr=(wsb+wsi)/2                             % m
% lr must be greater or equal than dr
dr= d3;
lr= dr;
Ar=lr*dr;                                            % m2
if lr==dr
    Le1=(ns*ns)*miu0*Lseg/32;                      % H
    Le1x=(nsx*nsx)*miu0*Lseg/32;                      % H
else
    Le1=(ns*ns)*(miu0*Lseg/(Ar^2))*(dr^4/32+ dr^3*(lr-dr)/16+ dr^2*(lr-dr)^2/64- dr*(lr-dr)^3/64+log(1+2*dr/(lr-dr))*(lr-dr)^4/128);
    Le1x=(nsx*nsx)*(miu0*Lseg/(Ar^2))*(dr^4/32+ dr^3*(lr-dr)/16+ dr^2*(lr-dr)^2/64- dr*(lr-dr)^3/64+log(1+2*dr/(lr-dr))*(lr-dr)^4/128);
end
% Pe2:
lew=lr; leo=0.1*lr; x=0.1*lr;                                             % m
y=0.2*dw; db=0.1*dw;
% Reluctance:
R= (leo+lew+x)/(db*Lseg*miu0)+(dw+db/2+y)/(2*x*Lseg*miu0)+(leo+lew+x)/(2*y*Lseg*miu0);          
Pel2=1/R;
Le2= (ns*ns*Pel2);                                  % H
Le2x= (nsx*nsx*Pel2);                               % H
% End Turn Inductance:
% End Leakage Inductance:
Le= (ns^2*miu0*Tc/8)*log(Tc^2*pi/(4*As));           % H
Lex= (nsx^2*miu0*Tc/8)*log(Tc^2*pi/(4*As));         % H
%Le= 2*Le1+2*Le2;                                   % H
GEN.ELECTRICAL.Le=Nsp*Le;                           % H
GEN.ELECTRICAL.Lex=Nsp*Lex;                         % H
% Reactances 
Xg= Nsp*we*Lg;                      % Ohms
Xgx= Nsp*we*Lgx;                    % Ohms
GEN.ELECTRICAL.Xg=Xg; 
GEN.ELECTRICAL.Xgx=Xgx; 
Xle1=Nsp*we*Le;                     % Ohms
Xle1x=Nsp*we*Lex;                   % Ohms
GEN.ELECTRICAL.Xe=Xle1;
GEN.ELECTRICAL.Xex=Xle1x;
Xls1=Nsp*we*Ls1;                    % Ohms
Xls1x=Nsp*we*Ls1x;                  % Ohms
GEN.ELECTRICAL.Xs=Xls1;
GEN.ELECTRICAL.Xs=Xls1x;
Xls1=Nsp*we*Ls1;                    % Ohms
Xls1x=Nsp*we*Ls1x;                  % Ohms
GEN.ELECTRICAL.Xs=Xls1;
GEN.ELECTRICAL.Xsx=Xls1x;

% Phase Inductance:
Lph= Nsp*(Lg+Ls1+ Le);              % H
Lphx= Nsp*(Lgx+Ls1x+ Lex);          % H
GEN.ELECTRICAL.Lph=Lph;             % H
GEN.ELECTRICAL.Lphx=Lphx;           % H
Xph= we*Lph;                        % Ohms
Xphx= we*Lphx;                      % Ohms
GEN.ELECTRICAL.Xph=Xph;             % H
GEN.ELECTRICAL.Xphx=Xphx;           % H
% STEADY STATE PARAMETERS:
Lmd=Nsp*Lg;
Lmdx=Nsp*Lgx;
Xmd=we*Lmd;
Xmdx=we*Lmdx;
Lmq=Nsp*Lg;
Lmqx=Nsp*Lgx;
Xmq=we*Lmq;
Xmqx=we*Lmqx;
Lq=Lg+Ls1;                          % H
Lqx=Lgx+Ls1x;                       % H
GEN.ELECTRICAL.Lq=Nsp*Lq;
GEN.ELECTRICAL.Lqx=Nsp*Lqx;  
Ld=Lg+Ls1;                          % H
Ldx=Lgx+Ls1x;                       % H
GEN.ELECTRICAL.Ld=Nsp*Ld; 
GEN.ELECTRICAL.Ldx=Nsp*Ldx; 
L0=Ls1;                             % H
L0x=Ls1x;                           % H
GEN.ELECTRICAL.L0=Nsp*L0;
GEN.ELECTRICAL.L0x=Nsp*L0x;
X0=Xls1;                            % Ohms
X0x=Xls1x;                          % Ohms
GEN.ELECTRICAL.X0=X0;
GEN.ELECTRICAL.X0x=X0x;
Xq=Xls1+Xmq;
Xqx=Xls1x+Xmqx;
GEN.ELECTRICAL.Xq=Xq;               % Ohms
GEN.ELECTRICAL.Xqx=Xqx;             % Ohms
Xd=Xls1+Xmd;
Xdx=Xls1x+Xmdx;
GEN.ELECTRICAL.Xd=Xd;               % Ohms
GEN.ELECTRICAL.Xdx=Xdx;             % Ohms
% Reluctancia del entrehierro (A-vuelta/Wb):
% Number of coils:
Ns1= ns*Ns/(2*pi);
Ns1x= ns*Ns/(2*pi);
Br=GEN.MAGNET.Br;
Brx=GEN.MAGNET.Brx;
Rg= (Rri/(miu0))*log(1+gef/(Rri+lm));   % A-vuelta/Wb
% Magnet Reluctance (Ri):
Rsp= (Rri/(miu0))*log(1+lm/Rri);
% Magnet Reluctance:
Rpm= (Rri/(miu0*miuR))*log(1+lm/Rri);
% Inductance Lqm (H):
% Lmq1= (6*Rri*L*(Ns1.^2)/(Nm.^2))*[(pi*(1-am)+sin(pi*am))/(Rsp+Rg)+(pi*am-sin(pi*am))/(Rg+Rpm)]
% % Reactance Xqm (Ohms):
% Xmq1=we*Lmq1                                                           % Ohms
% % Inductance Ldm:
% Lmd1= (6*Rri*L*(Ns1.^2)/(Nm.^2))*[(pi*(1-am)-sin(pi*am))/(Rsp+Rg)+(pi*am+sin(pi*am))/(Rg+Rpm)];
% % Reactance Xdm (Ohms):
% Xmd1=we*Lmd1; 
% lamdamp;
lamdam= (8*L*Rri*Ns1/(Nm*(Rpm+Rg))*(lm/(miu0*miuR)*Br*sin(pi*am/2)));
GEN.ELECTRICAL.lamdam=lamdam;                                           % Volts-s
%--------------------------------------------------------------------------
% Copper Weigth:
%--------------------------------------------------------------------------
%Ac=kcp*As;                                         % Conductor Area (m2)
GEN.WINDING.Ac=Ac;                                  % H
%dcu=sqrt(4*Ac/(pi));                               % m
GEN.WINDING.dcu=dcu;                                % m
Wcu= Ns*Ac*(L+pi*Tc)*Dcu;                           % kg
GEN.WINDING.Wcu=Wcu;                                % kg
% Copper volume:
Vcu=Ns*Ac*(L+pi*Tc);
GEN.WINDING.Vcu=Vcu;                                % kg
% Stator Volume:
%Vst=(pi*(Rso^2-Rsi^2)-Ns*As)*L*Kst;                % m3
Vst=(pi*(Rso^2-Rsi^2)-Ns*As)*L*Kst;                 % m3
GEN.STEEL.Vst=Vst;                                  % m3
Wst=Vst*rhobi;                                      % kg
GEN.STEEL.Wst=Wst;                                  % m3

% Rotor Volume:
%Vro=(pi*((Rro-lm)^2-Rri^2))*L*Kst;                % m3
Vro=(pi*((Rro-lm)^2-Rri^2))*L;                     % m3
GEN.STEEL.Vro=Vro;                                 % m3
Wro=Vro*rhobi;                                     % kg
GEN.STEEL.Wro=Wro;                                 % m3

% Total Weigth:
W=Wro+Wst+Wcu+Wm;                                  % kg
GEN.W.Data=W;                                      % kg

%--------------------------------------------------------------------------
% Estimated Inertia:
%--------------------------------------------------------------------------
Jrot=(1/2)*Wro*(Rri^2+(Rro-lm)^2)+(1/2)*Wm*(Rro^2+(Rro-lm)^2);                 % kg-m2
GEN.J.Data=Jrot;    
% else 
% GEN.ELECTRICAL.ns=ns;
% ns=ns*2;
% % Emax:
% Emax=cte*wm*ns;
% GEN.ELECTRICAL.Emax=Emax;
% % Peak Slot Current:
% Tin=GEN.Tem.Data;
% Is= Tin/(cte);
% GEN.ELECTRICAL.Is=Is;
% % Conductor current (A):
% Ic=Is/ns;
% GEN.ELECTRICAL.Ic=Ic;
% % Conductor Area (mm2):
% %Ac=Is/(Jmax*1e6);
% Ac=Is/Jmax/1e6;                      % m2
% %GEN.ELECTRICAL.Ac=Ac/1e6;           % m2
% GEN.ELECTRICAL.Ac=Ac;           % m2
% if Ac<As
%     kcp=Ac/As;
%     GEN.FACTORS.kcp=kcp;
% else
%     error('Aumentar los diámetros Dso, Dsi, Dri, Dro del PMSG');
% end
% % One conductor Area:
% Ac1=pi*dcu^2/4;                     % mm2
% GEN.ELECTRICAL.Ac1=Ac1;             % Ac1
% % Current allowed per single conductor:
% Ic1= Jmax*Ac1;                      % Amp
% GEN.ELECTRICAL.Ic1=Ic1;             % Amp
% % Number of parallel conductors:
% np= Ac*1e6/(ns*Ac1);                %
% GEN.ELECTRICAL.np=np;               % Amp
% % Phase current:
% Iph=Is/(Nph*ns);
% GEN.ELECTRICAL.Iph=Iph;
% % Peak conductor current density:
% Jc=Is/(kcp*As);                     % A/m2
% GEN.ELECTRICAL.Jc=Jc;
% % Peak slot flux density:
% miu0=4*pi*1e-7;     
% Bsmax= miu0*Is/ws;                       % H/m
% GEN.ELECTRICAL.Bsmax=Bsmax;              % Wb/m2
% % Slot Resistance (Ohms):
% % Conductor area (m2):
% %Ac=kcp*As;
% % Conductor heigth (m):
% wc=(wsi+wsb)/2;
% h= kcp*d3;                                  % m
% % Conductor diameter (m):
% %dcu=sqrt(4*Ac/(ns*pi));                    % m
% % Updating rhocu:
% rhocux=rhocu*(1+Beta*(Tx-20));
% GEN.WINDING.rhocux=rhocux;                  % Ohms/m
% % Skin depth:
% gam= sqrt(2*rhocux/(we*miu0));               % Copper profundity
% % DC Resistance:
% Acs=Ac/np;                                  % m2
% Rsx=rhocux*(ns*ns)*L/(Acs*np);              % Ohms
% Rs=rhocux*(ns*ns)*L/(Acs*np);               % Ohms
% % AC Resistance:
% % Estas fórmulas no funcionan bien:
% %Rac1= (ds*rhocu*L*dcu^2*ns^2/(9*gam^4*wsb))
% %deltae=(1/9)*(h^2/gam^2)*(dcu^2/gam^2)      %
% Racx= rhocux*(ns)*L/(pi*dcu*gam*np);         % Ohms
% Rac= rhocu*(ns)*L/(pi*dcu*gam*np);           % Ohms
% % Total Resistance:
% Rssx=sqrt(Rsx^2+Racx^2);                      % Ohms
% Rss=sqrt(Rsx^2+Racx^2);                       % Ohms
% %Rss1=sqrt(Rs^2+Rac1^2)                       % Ohms
% GEN.ELECTRICAL.Rs=Rssx;                       % Ohms
% % End-turn Resistance (Ohms):
% % DC:
% Rex=rhocux*(ns*ns)*pi*Tc/(2*Acs*np);         % Ohms
% Re=rhocu*(ns*ns)*pi*Tc/(2*Acs*np);           % Ohms
% % AC:
% Re_acx= rhocux*ns*pi*Tc/(pi*dcu*gam*np);     % Ohms
% Re_ac= rhocux*ns*pi*Tc/(pi*dcu*gam*np);      % Ohms
% % Total resistance:
% Reex=sqrt(Rex^2+Re_acx^2);                   % Ohms
% Ree=sqrt(Re^2+Re_ac^2);                      % Ohms
% GEN.ELECTRICAL.Re=Reex;                      % Ohms
% % Phase Resistance (Ohms):
% Rphx=Nsp*(Rssx+Reex);                        % Ohms
% GEN.ELECTRICAL.Rphx=Rphx;                    % Ohms
% GEN.ELECTRICAL.Rph=Rph;                      % Ohms
% % Phase Resistance (75°C):
% GEN.ELECTRICAL.Re=Ree;                      % Ohms
% % Resistance at 75°C:
% rhocu75=GEN.WINDING.rhocu75;                % Ohms-m
% Rph75=Rph*rhocu75/rhocu;
% GEN.ELECTRICAL.Rph75=Rph75;                 % Ohms
% %--------------------------------------------------------------------------
% % Inductance:
% %--------------------------------------------------------------------------
% % Air-Gap Inductance:
% kc=GEN.MAGNET.kc;
% gc=g+lm/miuR;
% GEN.MAGNET.gc=gc;
% gef= kc*gc;
% GEN.MAGNET.gef=gef;
% % Air-Gap Inductance (Inductance):
% Lg= (ns*ns)*miur*miu0*L*Tc*kd/(4*(lm+miur*kc*g));   % H
% GEN.ELECTRICAL.Lg=Nsp*Lg;  
% % Slot Leakage Inductance:
% % Ls= (ns*ns)*[miu0*d3^2*L/(3*(As))+ 2*miu0*d2*L/(ws+wsi) + miu0*d1*L/ws];  % H
% % Xls= Nsp*we*Ls;
% % Paths 5, 6 and 7:
% dw=0.765*d3;
% P1=miu0*gef*L/(wt+ws);
% P2=miu0*d1*L/(ws);
% P3=2*miu0*d2*L/(ws+wsi);
% P4=miu0*(d3-dw)*L/(wsi);
% P44=miu0*dw*L/(3*(wsi+wsb));
% P5=miu0*(d3/2)*log(1+pi*wt/wsi)/(3*pi);
% P6=miu0*log(1+pi*wt/wsi)*(d3/2)/pi;
% P7=miu0*log(1+pi*wt/wsi)*(d1+d2)/pi;
% Ls1=(ns*ns)*[P1+P2+P3+P4+P44+2*P5+2*P6+2*P7];       % H
% GEN.ELECTRICAL.Ls=Nsp*Ls1;  
% % End Leakage Inductance:
% % Segment length at the end of the turn:
% Lseg=2*pi*(Rsi+ds-0.5*dw)/Ns;                       % m
% % dr=d3; lr=(wsb+wsi)/2                             % m
% % lr must be greater or equal than dr
% dr= d3;
% lr= dr;
% Ar=lr*dr;                                            % m2
% if lr==dr
%     Le1=(ns*ns)*miu0*Lseg/32;                      % H
% else
%     Le1=(ns*ns)*(miu0*Lseg/(Ar^2))*(dr^4/32+ dr^3*(lr-dr)/16+ dr^2*(lr-dr)^2/64- dr*(lr-dr)^3/64+log(1+2*dr/(lr-dr))*(lr-dr)^4/128);
% end
% % Pe2:
% lew=lr; leo=0.1*lr; x=0.1*lr;                                             % m
% y=0.2*dw; db=0.1*dw;
% % Reluctance:
% R= (leo+lew+x)/(db*Lseg*miu0)+(dw+db/2+y)/(2*x*Lseg*miu0)+(leo+lew+x)/(2*y*Lseg*miu0);          
% Pel2=1/R;
% Le2= (ns*ns*Pel2);                                  % H
% % End Turn Inductance:
% % End Leakage Inductance:
% Le= (ns^2*miu0*Tc/8)*log(Tc^2*pi/(4*As));           % H
% %Le= 2*Le1+2*Le2;                                   % H
% GEN.ELECTRICAL.Le=Nsp*Le;                           % H
% 
% % Reactances 
% Xg= Nsp*we*Lg;                      % Ohms
% GEN.ELECTRICAL.Xg=Xg; 
% Xle1=Nsp*we*Le;                     % Ohms
% GEN.ELECTRICAL.Xe=Xle1;
% Xls1=Nsp*we*Ls1;                    % Ohms
% GEN.ELECTRICAL.Xs=Xls1;
% % Phase Inductance:
% Lph= Nsp*(Lg+Ls1+ Le);              % H
% GEN.ELECTRICAL.Lph=Lph;             % H
% Xph= we*Lph;                        % Ohms
% GEN.ELECTRICAL.Xph=Xph;             % H
% 
% % STEADY STATE PARAMETERS:
% Lmd=Nsp*Lg;
% Xmd=we*Lmd;
% Lmq=Nsp*Lg;
% Xmq=we*Lmq;
% Lq=Lg+Ls1;                          % H
% GEN.ELECTRICAL.Lq=Nsp*Lq;  
% Ld=Lg+Ls1;                          % H
% GEN.ELECTRICAL.Ld=Nsp*Ld; 
% L0=Ls1;                             % H
% GEN.ELECTRICAL.L0=Nsp*L0;
% X0=Xls1;                            % Ohms
% GEN.ELECTRICAL.X0=X0;
% Xq=Xls1+Xmq;
% GEN.ELECTRICAL.Xq=Xq;               % Ohms
% Xd=Xls1+Xmd;
% GEN.ELECTRICAL.Xd=Xd;               % Ohms
% % Reluctancia del entrehierro (A-vuelta/Wb):
% % Number of coils:
% Ns1= ns*Ns/(2*pi);
% Br=GEN.MAGNET.Br;
% Rg= (Rri/(miu0))*log(1+gef/(Rri+lm));   % A-vuelta/Wb
% % Magnet Reluctance (Ri):
% Rsp= (Rri/(miu0))*log(1+lm/Rri);
% % Magnet Reluctance:
% Rpm= (Rri/(miu0*miuR))*log(1+lm/Rri);
% % Inductance Lqm (H):
% % Lmq1= (6*Rri*L*(Ns1.^2)/(Nm.^2))*[(pi*(1-am)+sin(pi*am))/(Rsp+Rg)+(pi*am-sin(pi*am))/(Rg+Rpm)]
% % % Reactance Xqm (Ohms):
% % Xmq1=we*Lmq1                                                           % Ohms
% % % Inductance Ldm:
% % Lmd1= (6*Rri*L*(Ns1.^2)/(Nm.^2))*[(pi*(1-am)-sin(pi*am))/(Rsp+Rg)+(pi*am+sin(pi*am))/(Rg+Rpm)];
% % % Reactance Xdm (Ohms):
% % Xmd1=we*Lmd1; 
% % lamdamp;
% lamdam= (8*L*Rri*Ns1/(Nm*(Rpm+Rg))*(lm/(miu0*miuR)*Br*sin(pi*am/2)));
% GEN.ELECTRICAL.lamdam=lamdam;                                           % Volts-s
% %--------------------------------------------------------------------------
% % Copper Weigth:
% %--------------------------------------------------------------------------
% %Ac=kcp*As;                                         % Conductor Area (m2)
% GEN.WINDING.Ac=Ac;                                  % H
% %dcu=sqrt(4*Ac/(pi));                               % m
% GEN.WINDING.dcu=dcu;                                % m
% Wcu= Ns*Ac*(L+pi*Tc)*Dcu;                           % kg
% GEN.WINDING.Wcu=Wcu;                                % kg
% Vcu=Ns*Ac*(L+pi*Tc);
% GEN.WINDING.Vcu=Vcu;                                % m3
% % Stator Volume:
% %Vst=(pi*(Rso^2-Rsi^2)-Ns*As)*L*Kst;                 % m3
% Vst=(pi*(Rso^2-Rsi^2)-Ns*As)*L*Kst;                 % m3
% GEN.STEEL.Vst=Vst;                                  % m3
% Wst=Vst*rhobi;                                      % kg
% GEN.STEEL.Wst=Wst;                                  % m3
% 
% % Rotor Volume:
% %Vro=(pi*((Rro-lm)^2-Rri^2))*L*Kst;                 % m3
% Vro=(pi*((Rro-lm)^2-Rri^2))*L;                     % m3
% GEN.STEEL.Vro=Vro;                                 % m3
% Wro=Vro*rhobi;                                     % kg
% GEN.STEEL.Wro=Wro;                                 % m3
% 
% % Total Weigth:
% W=Wro+Wst+Wcu+Wm;                                  % kg
% GEN.W.Data=W;                                      % kg
% 
% %--------------------------------------------------------------------------
% % Estimated Inertia:
% %--------------------------------------------------------------------------
% Jrot=(1/2)*Wro*(Rri^2+(Rro-lm)^2)+(1/2)*Wm*(Rro^2+(Rro-lm)^2);                 % kg-m2
% GEN.J.Data=Jrot;  
% end

