% ------------------------------------------------------------------------
% Main script for computing SDOS 
% ------------------------------------------------------------------------
clear all
close all
clc

% % Parameters
a = 1;                                % lattice constant
c = 3e8;                              % speed of light
nf = 21;                              % number of frequency
nk = 21;                              % number of wavevector
freq = c/a*linspace(0.5, 0.7, nf);    % normalized frequency
omega = 2*pi*freq;                    % angular frequency
k = linspace(-pi, pi, nk)/a;          % wave vector
phase = exp(1i*k*a);                  % bloch phase
SDOS = zeros(nf, nk);                 % surface density of states
eta = 1e-3;                           % imaginary part of frequency
posTol = a*1e-3;                      % position tolerance

% Load FEM information and require COMSOL-MATLAB environment
[mtrE, mtrK, no2xyz] = loadFEM;

% Find specific nodes indexes
noIdx_cell0=find(no2xyz(2,:)<=1*a-posTol & no2xyz(2,:)>=-0.5*a-posTol);
noIdx_cell1=find(no2xyz(2,:)<=-0.5*a-posTol & no2xyz(2,:)>=-1.5*a-posTol);
noIdx_cell2=find(no2xyz(2,:)<=-1.5*a-posTol & no2xyz(2,:)>=-2.5*a-posTol);
noIdx_cell0TopEdge=find(abs(no2xyz(2,:)-(1*a))<posTol);

% Reconstruct node indexes
[~,noIdx_tmpCell]=min(pdist2(no2xyz(:,noIdx_cell1)',no2xyz(:,noIdx_cell2)'+[0,1*a])');
noIdx_cell2=noIdx_cell2(:,noIdx_tmpCell);

% Reconstruct node positions
noIdx_tmpAll=[noIdx_cell0TopEdge,noIdx_cell0,noIdx_cell1,noIdx_cell2];
no2xyz=no2xyz(:,noIdx_tmpAll);

% Reconstruct mass and stiffness matrices
mtrE=mtrE(noIdx_tmpAll(1:end),noIdx_tmpAll(1:end));
mtrK=mtrK(noIdx_tmpAll(1:end),noIdx_tmpAll(1:end));

% Node indexes for cell0 top edge (PEC)
noIdx_cell0TopEdge = 1:length(noIdx_cell0TopEdge);

% Node indexes outside cell0 top edge (PEC)
noIdx_int = length(noIdx_cell0TopEdge)+1:length(no2xyz);

% Node indexes for left and right edges (PBC)
noIdx_leftEdge=find(abs(no2xyz(1,:)-(-0.5*a))<posTol & no2xyz(2,:)<=1*a-posTol);
noIdx_rightEdge=find(abs(no2xyz(1,:)-(2.5*a))<posTol & no2xyz(2,:)<=1*a-posTol);
[~,noIdx_tmpEdge]=min(pdist2(no2xyz(:,noIdx_leftEdge)',no2xyz(:,noIdx_rightEdge)'-[3*a,0])');
noIdx_rightEdge=noIdx_rightEdge(:,noIdx_tmpEdge);

% Define DOFs for truncating matrices
dof_cell0LeftEdge=length(find(abs(no2xyz(1,:)-(-0.5*a))<posTol & no2xyz(2,:)<=1*a-posTol & no2xyz(2,:)>=-0.5*a-posTol));
dof_cell1LeftEdge=length(find(abs(no2xyz(1,:)-(-0.5*a))<posTol & no2xyz(2,:)<=-0.5*a-posTol & no2xyz(2,:)>=-1.5*a-posTol));
dof_cell2LeftEdge=dof_cell1LeftEdge;
dof_cell0Int=length(noIdx_cell0)-dof_cell0LeftEdge;
dof_cell1Int=length(noIdx_cell1)-dof_cell1LeftEdge;
dof_cell2Int=length(noIdx_cell2)-dof_cell2LeftEdge;

dof_all = length(no2xyz);
dof_cell0TopEdge = length(noIdx_cell0TopEdge); 
dof_leftEdge = length(noIdx_leftEdge);

% Compute SDOS 
for mf = 1:nf
    for mk = 1:nk
        
        % Initialize the matrix for imposing boundary conditions
        Pn=eye(dof_all);
        
        % Impose PEC boundary condition
        for i = 1:dof_cell0TopEdge
            Pn(noIdx_cell0TopEdge(i), noIdx_cell0TopEdge(i)) = 0;
        end
        
        % Impose periodic boundary conditions
        for i = 1:dof_leftEdge
            Pn(noIdx_leftEdge(i), noIdx_leftEdge(i)) = 0;
            Pn(noIdx_leftEdge(i), noIdx_rightEdge(i)) = phase(mk);
        end
        
        % Reduce matrix dimension
        Pn = Pn(any(Pn, 2), any(Pn, 1));
        Pn = sparse(Pn);
        
        % Construct eliminated mass and stiffness matrices
        mtrKc=Pn'*mtrK(noIdx_int,noIdx_int)*Pn;
        mtrEc=Pn'*mtrE(noIdx_int,noIdx_int)*Pn;
        
        Z = (1-1i*eta)^2 * omega(mf)^2 * mtrEc - mtrKc;
        Z00=Z(1:dof_cell0Int,1:dof_cell0Int);
        Z01=Z(1:dof_cell0Int,dof_cell0Int+1:dof_cell0Int+dof_cell1Int);
        Z10=Z(dof_cell0Int+1:dof_cell0Int+dof_cell1Int,1:dof_cell0Int);
        Z11=Z(dof_cell0Int+1:dof_cell0Int+dof_cell1Int,dof_cell0Int+1:dof_cell0Int+dof_cell1Int);
        Z12=Z(dof_cell0Int+1:dof_cell0Int+dof_cell1Int,dof_cell0Int+dof_cell1Int+1:dof_cell0Int+dof_cell1Int+dof_cell2Int);
        Z21=Z(dof_cell0Int+dof_cell1Int+1:dof_cell0Int+dof_cell1Int+dof_cell2Int,dof_cell0Int+1:dof_cell0Int+dof_cell1Int);
        
        % Implement CRM or TMM
        G00=crmSGF(Z11,Z12,Z21,Z00,Z01,Z10);
        % G00=tmmSGF(Z11,Z12,Z21,Z00,Z01,Z10);
        SDOS(mf,mk) = omega(mf) * imag(trace(G00));
    end
    disp(mf);
end