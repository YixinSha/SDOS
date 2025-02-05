% ------------------------------------------------------------------------
% Calculate G00 using cyclic reduction method 
% ------------------------------------------------------------------------
function [G00] = crmSGF(Z11pos,Z12pos,Z21pos,Z11neg,Z12neg,Z21neg,Z00,Z01pos,Z10pos,Z01neg,Z10neg)

maxSteps = 20;  % maximum iterations
convTol = 1e-2; % convergence tolerance
divTol = 1e3;   % divergence tolerance

AlphaPos(:,:,1)=full(Z12pos);
AlphaPosS(:,:,1)=full(Z01pos);
BetaPos(:,:,1)=full(Z21pos);
BetaPosS(:,:,1)=full(Z10pos);
AlphaNeg(:,:,1)=full(Z12neg);
AlphaNegS(:,:,1)=full(Z01neg);
BetaNeg(:,:,1)=full(Z21neg);
BetaNegS(:,:,1)=full(Z10neg);
ZetaPos(:,:,1)=full(Z11pos);
ZetaNeg(:,:,1)=full(Z11neg);
ZetaS(:,:,1)=full(Z00);

for i=1:maxSteps
    
    ZetaPos_Div_AlphaPos=ZetaPos(:,:,1)\AlphaPos(:,:,1);
    ZetaPos_Div_BetaPos=ZetaPos(:,:,1)\BetaPos(:,:,1);
    ZetaPos_Div_BetaPosS=ZetaPos(:,:,1)\BetaPosS(:,:,1);
    
    ZetaNeg_Div_AlphaNeg=ZetaNeg(:,:,1)\AlphaNeg(:,:,1);
    ZetaNeg_Div_BetaNeg=ZetaNeg(:,:,1)\BetaNeg(:,:,1);
    ZetaNeg_Div_BetaNegS=ZetaNeg(:,:,1)\BetaNegS(:,:,1);
    
    AlphaPos(:,:,2)=AlphaPos(:,:,1)*(ZetaPos_Div_AlphaPos);
    AlphaPosS(:,:,2)=AlphaPosS(:,:,1)*(ZetaPos_Div_AlphaPos);
    BetaPos(:,:,2)=BetaPos(:,:,1)*(ZetaPos_Div_BetaPos);
    BetaPosS(:,:,2)=BetaPos(:,:,1)*(ZetaPos_Div_BetaPosS);
    
    AlphaNeg(:,:,2)=AlphaNeg(:,:,1)*(ZetaNeg_Div_AlphaNeg);
    AlphaNegS(:,:,2)=AlphaNegS(:,:,1)*(ZetaNeg_Div_AlphaNeg);
    BetaNeg(:,:,2)=BetaNeg(:,:,1)*(ZetaNeg_Div_BetaNeg);
    BetaNegS(:,:,2)=BetaNeg(:,:,1)*(ZetaNeg_Div_BetaNegS);
    
    ZetaPos(:,:,2)=ZetaPos(:,:,1)-AlphaPos(:,:,1)*(ZetaPos_Div_BetaPos)-BetaPos(:,:,1)*(ZetaPos_Div_AlphaPos);
    ZetaNeg(:,:,2)=ZetaNeg(:,:,1)-AlphaNeg(:,:,1)*(ZetaNeg_Div_BetaNeg)-BetaNeg(:,:,1)*(ZetaNeg_Div_AlphaNeg);
    ZetaS(:,:,2)=ZetaS(:,:,1)-AlphaPosS(:,:,1)*(ZetaPos_Div_BetaPosS)-AlphaNegS(:,:,1)*(ZetaNeg_Div_BetaNegS);
    
    dZetaS = ZetaS(:, :, 2) - ZetaS(:, :, 1);
    error = norm(dZetaS, "fro") / norm(ZetaS(:, :, 1), "fro");
    
    AlphaPos(:,:,1)=AlphaPos(:,:,2);
    AlphaPosS(:,:,1)=AlphaPosS(:,:,2);
    BetaPos(:,:,1)=BetaPos(:,:,2);
    BetaPosS(:,:,1)=BetaPosS(:,:,2);
    AlphaNeg(:,:,1)=AlphaNeg(:,:,2);
    AlphaNegS(:,:,1)=AlphaNegS(:,:,2);
    BetaNeg(:,:,1)=BetaNeg(:,:,2);
    BetaNegS(:,:,1)=BetaNegS(:,:,2);
    ZetaPos(:,:,1)=ZetaPos(:,:,2);
    ZetaNeg(:,:,1)=ZetaNeg(:,:,2);
    ZetaS(:,:,1)=ZetaS(:,:,2);
    
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