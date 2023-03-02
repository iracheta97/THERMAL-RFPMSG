function[GEN]= TOPOLOGICAL_CONSTRAINT(PMSG)

% Inputs:
% PMSG: Structure Data of the PMSG

% Outputs:
% GEN: Updating the Structure Data of the PMSG
GEN=PMSG;
% Electrical Frequency (Hz):
fe= PMSG.fe.Data;

% Rated Speed (rpm):
Sr= PMSG.Sr.Data;           % rpm

%--------------------------------------------------------------------------
% Number of poles:
%--------------------------------------------------------------------------
Nm=GEN.CONSTRAINTS.Nm; 

% Nm must be even and integer
% Verify if Nm is integer:
X= Nm-floor(Nm);

if X==0     % Integer
    % Vefify if Nm is an even number
    Y= rem(Nm,2);
    if Y==0
    % Number of poles:
    GEN.CONSTRAINTS.Nm=Nm;
    else
    % Add one poles
    Nm=Nm+1;
    % Number of poles:
    GEN.CONSTRAINTS.Nm=Nm;
    % Recalculate Sr:
    Sr= 120*fe/Nm;
    GEN.Sr.Data= Sr;           % rpm
    end
else      % Not Integer:
    % Round to the floor integer
    Nm=floor(Nm);
    % Verify if Nm is an even number
    Y=rem(Nm,2);
    if Y==0
    % Number of poles:
    GEN.CONSTRAINTS.Nm=Nm;
    % Recalculate Sr:
    Sr= 120*fe/Nm;
    GEN.Sr.Data= Sr;           % rpm
    else
    % Add one poles
    Nm=Nm+1;
    % Number of poles:
    GEN.CONSTRAINTS.Nm=Nm;
    % Recalculate Sr:
    Sr= 120*fe/Nm;
    GEN.Sr.Data= Sr;           % rpm   
    end
end
%--------------------------------------------------------------------------
% Number of Slots:
%--------------------------------------------------------------------------
% Number of phases:
Nph=GEN.CONSTRAINTS.Nph;
% Number of poles:
Nm=GEN.CONSTRAINTS.Nm;
% Number of slots:
Ns=GEN.CONSTRAINTS.Ns;                     % Multiple Integer of Nph 
GEN.CONSTRAINTS.Ns=Ns;
% Number of slots per phase (Nsp>Nm):
%Nsp= Nm+2;                      % where Nsp is an even integer number of slots per phase
Nsp=Ns/Nph;
GEN.CONSTRAINTS.Nsp=Nsp;
% Number of slots per pole per phase:
Nspp= Ns/(Nm*Nph);              % Nspp >=1 (This variable can be fractional)
GEN.CONSTRAINTS.Nspp=Nspp;
% Number of slots per magnet>
Nsm= Nspp*Nph;                  % or Ns/Nm (Nsm can be fractional)
GEN.CONSTRAINTS.Nsm=Nsm;
%



