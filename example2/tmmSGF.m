% ------------------------------------------------------------------------
% Calculate G00 using transfer matrix method 
% ------------------------------------------------------------------------
function [G00] = tmmSGF(Z11, Z12, Z21, Z00, Z01, Z10)
% Arguments:                                                  
%    Z00 = intra-coupling matrix
%    Z01 = inter-coupling matrix
%    Z10 = inter-coupling matrix
%    Z11 = intra-coupling matrix
%    Z12 = inter-coupling matrix
%    Z21 = inter-coupling matrix
% Returns:                                                
%    G00 = surface Green's function
sizeZ11 = length(Z11); % size of block matrix

% Construct T1 and T2 matrices
T1 = [zeros(sizeZ11), eye(sizeZ11); -Z12, zeros(sizeZ11)];
T2 = [eye(sizeZ11), zeros(sizeZ11); Z11, Z21];

% Do generalized eigenanalysis
[eigVtr, eigVal] = eig(full(T2), full(T1));

% Order eigenvalues and eigenvectors
[~, eigIdx] = sort(abs(diag(eigVal)));
eigVtr = eigVtr(:, eigIdx);
S1 = eigVtr(sizeZ11+1:end, 1:sizeZ11);
S2 = eigVtr(1:sizeZ11, 1:sizeZ11);


% Calculate surface Green's function
G00 = inv(Z00 - Z01 / (Z11 + Z12 * S2 / S1) * Z10);

end