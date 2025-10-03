function [ lambda,minVal,lambda_range,creso_vals ] = cresoJD( matU, vectS, vectB )
%CRESOJD Perform CRESO calculation of regulation parameter for Tikhonov
% Inputs:
%   matU: SVD U matrix of A matrix
%   vectS: Singular values
%   vectB: vector B
%
% JD 20/11/14 - see Johnston PR et al. IEEE Trans BioMed Eng 1997
    

%% Input verification
    if min(size(vectB)) > 1, error('Incorrect input B vector'); end

%% Precalculation and settings
    vectAlpha = matU'*vectB;
    vectAlpha = vectAlpha(:);
    S = vectS(:);
    numElem = length(S);
    oneVector = ones(numElem,1);
    
    % Set minimal and maximal regulation parameters
    minReg = 10^-8; maxReg = 10^-2;
     
%% Define CRESO function
    function [cresoval] = local_creso(t)
        t = t.*oneVector;
        cresoval = (S .* vectAlpha ./ (S .^ 2 + t)) .^ 2;
        cresoval = cresoval .* (oneVector - (4 * t ./(S .^ 2 + t)));
        cresoval = -sum(cresoval, 1);
    end

%% Find local minimum of CRESO function
   [lambda,minVal]  = fminbnd(@local_creso, minReg, maxReg);

    

%% Plot CRESO over log-spaced lambda range
lambda_range = logspace(log10(minReg), log10(maxReg), 500);
creso_vals = zeros(size(lambda_range));

for i = 1:length(lambda_range)
    creso_vals(i) = local_creso(lambda_range(i));
end


end

