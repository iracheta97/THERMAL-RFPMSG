function[GEN]= PMSG_GEOMETRY(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
% Topological Constraints:
Nph=GEN.CONSTRAINTS.Nph;
Nm=GEN.CONSTRAINTS.Nm;
Ns=GEN.CONSTRAINTS.Ns;
Nspp=GEN.CONSTRAINTS.Nspp;
% -------------------------------------------------------------------------
% OPTIMIZING MAGNET FRACTION: (FOR AS=0.5)
%--------------------------------------------------------------------------
% Optimal Magnet Fraction am= Tm/Tp when as=ws/Ts=0.5
% Decreases Tm and lm to avoid Cogging Torque 
for n=1:10
    am(n)=(n+0.14)/(Nspp*Nph);
     if am(n)>1
        break;
     end  
end
% Optimal value for as=0.5
am2=0.6711;
%am2=am(end-1);
GEN.GEOMETRY.am=am2;
% Flux Concentration Factor:
Cphi= 2*am2/(1+am2);
% Rated Speed (rpm):
GEN.GEOMETRY.Cphi=Cphi;
%--------------------------------------------------------------------------
% Geometric Relationships:
%--------------------------------------------------------------------------
% Outside Rotor Radius:
Rro=GEN.GEOMETRY.Rro;          % m
% Gap:
g= GEN.GEOMETRY.g;             % m
% Number of slots:
Ns=GEN.CONSTRAINTS.Ns;         % m
% Nspp:
Nspp=GEN.CONSTRAINTS.Nspp;     % m
% Inside Stator Radius:
Rsi= Rro + g;               % m
GEN.GEOMETRY.Rsi=Rsi;       % m
% Pole pitch (Tp):
Tp=2*pi*Rro/Nm;
GEN.GEOMETRY.Tp=Tp;         % m
% Angular Pole-Pitch in Mechanical Radians
Tetap=2*pi/Nm;                      % radians
GEN.GEOMETRY.Tetap=Tetap;           % radians
% Slot pitch (Ts):
Ts=2*pi*Rsi/Ns;
GEN.GEOMETRY.Ts=Ts;         % m
% Angular Slot-Pitch in Mechanical Radians:
Tetas=2*pi/Ns;                      % radians
GEN.GEOMETRY.Tetas=Tetas;           % radians
% Angular Slot-Pitch in Electrical Radians:
Tetase=pi/(Nspp*Nph);               % radians
GEN.GEOMETRY.Tetase=Tetase;         % radians
% Coil-pitch Fraction:
if Nspp>=1
        acp= floor(Nspp)/Nspp;
else
        acp= 1/Nspp;
end
GEN.GEOMETRY.acp=acp;              % radians
% Coil-Pitch:
Tc= acp*Tp;                       % m
GEN.GEOMETRY.Tc=Tc;               % m
% Coil-Pitch in radians:
Tetace= acp*pi;                    % radians
GEN.GEOMETRY.Tetace=Tetace;        % radians

