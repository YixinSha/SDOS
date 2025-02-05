% ------------------------------------------------------------------------
% Calculate G00 using transfer matrix method 
% ------------------------------------------------------------------------
function [G00] = tmmSGF(Z11pos,Z12pos,Z21pos,Z11neg,Z12neg,Z21neg,Z00,Z01pos,Z10pos,Z01neg,Z10neg)

sizeZ11pos = length(Z11pos);
sizeZ11neg = length(Z11neg);

T1pos=[zeros(sizeZ11pos),eye(sizeZ11pos);-Z12pos,zeros(sizeZ11pos)];
T2pos=[eye(sizeZ11pos),zeros(sizeZ11pos);Z11pos,Z21pos];

T1neg=[zeros(sizeZ11neg),eye(sizeZ11neg);-Z12neg,zeros(sizeZ11neg)];
T2neg=[eye(sizeZ11neg),zeros(sizeZ11neg);Z11neg,Z21neg];

% Do generalized eigenanalysis
[eigVtr, eigVal] = eig(full(T2pos), full(T1pos));

% Order eigenvalues and eigenvectors
[~, eigIdx] = sort(abs(diag(eigVal)));
eigVtr = eigVtr(:, eigIdx);
S1pos = eigVtr(sizeZ11pos+1:end, 1:sizeZ11pos);
S2pos = eigVtr(1:sizeZ11pos, 1:sizeZ11pos);

% Do generalized eigenanalysis
[eigVtr, eigVal] = eig(full(T2neg), full(T1neg));

% Order eigenvalues and eigenvectors
[~, eigIdx] = sort(abs(diag(eigVal)));
eigVtr = eigVtr(:, eigIdx);
S1neg = eigVtr(sizeZ11neg+1:end, 1:sizeZ11neg);
S2neg = eigVtr(1:sizeZ11neg, 1:sizeZ11neg);

G00=inv(Z00-Z01pos/(Z11pos+Z12pos*S2pos/S1pos)*Z10pos-Z01neg/(Z11neg+Z12neg*S2neg/S1neg)*Z10neg);