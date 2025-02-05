% ------------------------------------------------------------------------
% Calculate G00 using cyclic reduction method 
% ------------------------------------------------------------------------
function [G00] = crmSGF(Z11, Z12, Z21, Z00, Z01, Z10)
% Arguments:                                                  
%    Z00 = intra-coupling matrix
%    Z01 = inter-coupling matrix
%    Z10 = inter-coupling matrix
%    Z11 = intra-coupling matrix
%    Z12 = inter-coupling matrix
%    Z21 = inter-coupling matrix
% Returns:                                                
%    G00 = surface Green's function 

maxSteps = 20;  % maximum iterations
convTol = 1e-2; % convergence tolerance
divTol = 1e3;   % divergence tolerance

% Initialize viriables
Alpha(:, :, 1) = full(Z12);
AlphaS(:, :, 1) = full(Z01);
Beta(:, :, 1) = full(Z21);
BetaS(:, :, 1) = full(Z10);
Zeta(:, :, 1) = full(Z11);
ZetaS(:, :, 1) = full(Z00);

for i = 1:maxSteps
    
    % Update viriables
    ZetaDivAlpha = Zeta(:, :, 1) \ Alpha(:, :, 1);
    ZetaDivBeta = Zeta(:, :, 1) \ Beta(:, :, 1);
    ZetaDivBetaS = Zeta(:, :, 1) \ BetaS(:, :, 1);
    
    Alpha(:, :, 2) = Alpha(:, :, 1) * (ZetaDivAlpha);
    AlphaS(:, :, 2) = AlphaS(:, :, 1) * (ZetaDivAlpha);
    Beta(:, :, 2) = Beta(:, :, 1) * (ZetaDivBeta);
    BetaS(:, :, 2) = Beta(:, :, 1) * (ZetaDivBetaS);
    Zeta(:, :, 2) = Zeta(:, :, 1) - Alpha(:, :, 1) * (ZetaDivBeta) - Beta(:, :, 1) * (ZetaDivAlpha);
    ZetaS(:, :, 2) = ZetaS(:, :, 1) - AlphaS(:, :, 1) * (ZetaDivBetaS);
    
    dZetaS = ZetaS(:, :, 2) - ZetaS(:, :, 1);
    error = norm(dZetaS, "fro") / norm(ZetaS(:, :, 1), "fro");

    Alpha(:, :, 1) = Alpha(:, :, 2);
    AlphaS(:, :, 1) = AlphaS(:, :, 2);
    Beta(:, :, 1) = Beta(:, :, 2);
    BetaS(: ,:, 1) = BetaS(:, :, 2);
    Zeta(:, :, 1) = Zeta(:, :, 2);
    ZetaS(:, :, 1) = ZetaS(:, :, 2);

    % Check for convergence and divergence
    if error < convTol
        break;
    elseif error > divTol
        disp('error');
        break;
    else
        continue;
    end

end

% Calculate surface Green's function
G00 = inv(ZetaS(:, :, 1));

end