% ------------------------------------------------------------------------
% Calculate G00 using transfer matrix method 
% ------------------------------------------------------------------------
function [G00] = tmmSGF(Z00, Z01, Z10)
% Arguments:                                                  
%    Z00 = intra-coupling matrix
%    Z01 = inter-coupling matrix
%    Z10 = inter-coupling matrix
% Returns:                                                
%    G00 = surface Green's function                      

sizeZ00 = length(Z00); % size of block matrix

% Construct T1 and T2 matrices
T1 = [zeros(sizeZ00), eye(sizeZ00); -Z01, zeros(sizeZ00)];
T2 = [eye(sizeZ00), zeros(sizeZ00); Z00, Z10];

% Do generalized eigenanalysis
[eigVtr, eigVal] = eig(full(T2), full(T1));

% Order eigenvalues and eigenvectors
[~, eigIdx] = sort(abs(diag(eigVal)));
eigVtr = eigVtr(:, eigIdx);
S1 = eigVtr(sizeZ00+1:end, 1:sizeZ00);
S2 = eigVtr(1:sizeZ00, 1:sizeZ00);


% Calculate surface Green's function
G00 = inv(Z00 + Z01 * S2 / S1);

end