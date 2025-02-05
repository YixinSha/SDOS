% ------------------------------------------------------------------------
% Main script for computing SDOS 
% ------------------------------------------------------------------------
clear all
close all
clc

% Parameters
a = 1;                                % lattice constant
c = 3e8;                              % speed of light
nf = 21;                              % number of frequency
nk = 21;                              % number of wavevector
freq = c/a*linspace(0.15, 0.35, nf);  % normalized frequency
omega = 2*pi*freq;                    % angular frequency
k = linspace(-pi, pi, nk)/a;          % wave vector
phase = exp(1i*k*a);                  % bloch phase
SDOS = zeros(nf, nk);                 % surface density of states
eta = 1e-3;                           % imaginary part of frequency
posTol = a*1e-3;                      % position tolerance

% Load FEM information and require COMSOL-MATLAB environment
[mtrE, mtrK, no2xyz] = loadFEM;

% Find specific nodes indexes
noIdx_cellNeg2=find(no2xyz(2,:)>=-3*a/sqrt(3)/4+a*sqrt(3)/2+a*sqrt(3)/2+posTol);
noIdx_cellNeg1=find(no2xyz(2,:)<=-3*a/sqrt(3)/4+a*sqrt(3)/2+a*sqrt(3)/2+posTol & no2xyz(2,:)>=-3*a/sqrt(3)/4+a*sqrt(3)/2+posTol);
noIdx_cell0=find(no2xyz(2,:)<=-3*a/sqrt(3)/4+a*sqrt(3)/2+posTol & no2xyz(2,:)>=-3*a/sqrt(3)/4-posTol);
noIdx_cellPos1=find(no2xyz(2,:)<=-3*a/sqrt(3)/4-posTol & no2xyz(2,:)>=-3*a/sqrt(3)/4-a*sqrt(3)/2-posTol);
noIdx_cellPos2=find(no2xyz(2,:)<=-3*a/sqrt(3)/4-a*sqrt(3)/2-posTol);

% Reconstruct node indexes in negative direction
[~,noIdx_tmpCellNeg]=min(pdist2(no2xyz(:,noIdx_cellNeg2)',no2xyz(:,noIdx_cellNeg1)'+[-a/2,a*sqrt(3)/2])');
noIdx_cellNeg1=noIdx_cellNeg1(:,noIdx_tmpCellNeg);

% Reconstruct node indexes in positive direction
[~,noIdx_tmpCellPos]=min(pdist2(no2xyz(:,noIdx_cellPos1)',no2xyz(:,noIdx_cellPos2)'+[-a/2,a*sqrt(3)/2])');
noIdx_cellPos2=noIdx_cellPos2(:,noIdx_tmpCellPos);

% Reconstruct node positions
noIdx_tmpAll=[noIdx_cellNeg2,noIdx_cellNeg1,noIdx_cell0,noIdx_cellPos1,noIdx_cellPos2];
no2xyz=no2xyz(:,noIdx_tmpAll);

% Reconstruct mass and stiffness matrices
mtrE=mtrE(noIdx_tmpAll(1:end),noIdx_tmpAll(1:end));
mtrK=mtrK(noIdx_tmpAll(1:end),noIdx_tmpAll(1:end));

% Node indexes for left and right edges (PBC)
noIdx_leftEdge=find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol);
noIdx_rightEdge=find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)+sqrt(3)/2*a))<posTol);
[~,noIdx_tmpEdge]=min(pdist2(no2xyz(:,noIdx_leftEdge)',no2xyz(:,noIdx_rightEdge)'-[1*a,0])');
noIdx_rightEdge=noIdx_rightEdge(:,noIdx_tmpEdge);

% Define DOFs for truncating matrices
dof_cellNeg2LeftEdge=length(find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol & no2xyz(2,:)>=(-3*a/sqrt(3)/4+a*sqrt(3)/2+a*sqrt(3)/2)+posTol));
dof_cellNeg1LeftEdge=dof_cellNeg2LeftEdge;
dof_cell0LeftEdge=length(find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol & no2xyz(2,:)<=(-3*a/sqrt(3)/4+a*sqrt(3)/2)+posTol & no2xyz(2,:)>=(-3*a/sqrt(3)/4)-posTol));
dof_cellPos2LeftEdge=length(find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol & no2xyz(2,:)<=(-3*a/sqrt(3)/4-a*sqrt(3)/2)-posTol));
dof_cellPos1LeftEdge=dof_cellPos2LeftEdge;

dof_cellNeg2Int=length(noIdx_cellNeg2)-dof_cellNeg2LeftEdge;
dof_cellNeg1Int=length(noIdx_cellNeg1)-dof_cellNeg1LeftEdge;
dof_cell0Int=length(noIdx_cell0)-dof_cell0LeftEdge;
dof_cellPos1Int=length(noIdx_cellPos1)-dof_cellPos1LeftEdge;
dof_cellPos2Int=length(noIdx_cellPos2)-dof_cellPos2LeftEdge;

dof_all = length(no2xyz);
dof_leftEdge = length(noIdx_leftEdge);

% Compute SDOS 
for mf = 1:nf
    for mk = 1:nk

        % Initialize the matrix for imposing boundary conditions
        Pn=eye(dof_all);

        % Impose periodic boundary conditions
        for i=1:dof_leftEdge
            Pn(noIdx_leftEdge(i),noIdx_leftEdge(i))=0;
            Pn(noIdx_leftEdge(i),noIdx_rightEdge(i))=phase(mk);
        end

        % Reduce matrix dimension
        Pn = Pn(any(Pn, 2), any(Pn, 1));
        Pn = sparse(Pn);
        
        % Construct eliminated mass and stiffness matrices
        mtrKc=Pn'*mtrK*Pn;
        mtrEc=Pn'*mtrE*Pn;
        
        % Construct block Z matrix
        Z = (1-1i*eta)^2 * omega(mf)^2 * mtrEc - mtrKc;
        
        Z00=Z(dof_cellNeg2Int+dof_cellNeg1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int,dof_cellNeg2Int+dof_cellNeg1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int);
        Z01pos=Z(dof_cellNeg2Int+dof_cellNeg1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int,dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int);
        Z10pos=Z(dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int,dof_cellNeg2Int+dof_cellNeg1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int);
        Z11pos=Z(dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int,dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int);
        Z12pos=Z(dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int,dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int+dof_cellPos2Int);
        Z21pos=Z(dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int+dof_cellPos2Int,dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int+dof_cellPos1Int);
        Z01neg=Z(dof_cellNeg2Int+dof_cellNeg1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int,dof_cellNeg2Int+1:dof_cellNeg2Int+dof_cellNeg1Int);
        Z10neg=Z(dof_cellNeg2Int+1:dof_cellNeg2Int+dof_cellNeg1Int,dof_cellNeg2Int+dof_cellNeg1Int+1:dof_cellNeg2Int+dof_cellNeg1Int+dof_cell0Int);
        Z11neg=Z(dof_cellNeg2Int+1:dof_cellNeg2Int+dof_cellNeg1Int,dof_cellNeg2Int+1:dof_cellNeg2Int+dof_cellNeg1Int);
        Z12neg=Z(dof_cellNeg2Int+1:dof_cellNeg2Int+dof_cellNeg1Int,1:dof_cellNeg2Int);
        Z21neg=Z(1:dof_cellNeg2Int,dof_cellNeg2Int+1:dof_cellNeg2Int+dof_cellNeg1Int);
        
        % Implement CRM or TMM
        G00=crmSGF(Z11pos,Z12pos,Z21pos,Z11neg,Z12neg,Z21neg,Z00,Z01pos,Z10pos,Z01neg,Z10neg);
        % G00=tmmSGF(Z11pos,Z12pos,Z21pos,Z11neg,Z12neg,Z21neg,Z00,Z01pos,Z10pos,Z01neg,Z10neg);
        SDOS(mf,mk) = omega(mf) * imag(trace(G00));
    end
    disp(mf);
end