function y = robust_zscore(x1,x2)
% [] W.A. Stahel. Robuste Schatzungen: infinitesimale Optimalit¨at und Schatzungen von Kovarianzmatrizen. PhD thesis, ETH Zurich, 1981.
% [] D.L. Donoho. Breakdown properties of multivariate location estimators. Qualifying paper, Harvard University, Boston, 1982.

% x1, x2 are either vectirs or matrices
% If a matrix (eg Nsubj x Neeg), then oparates along the columns

if nargin == 1
    % disp('Robust Z-score using patient data.');
    y = (x1-median(x1)) ./ (1.4826*mad(x1,1));     % median absolute deviation
    % y = (x1-median(x1)) ./(1.253314*mad(x1,0));  % mean absolute deviation
else
    % disp('Robust Z-score using control data.');
    y = (x1-median(x2)) ./ (1.4826*mad(x2,1));      % median absolute deviation
    % y = (x1-median(x2)) ./ (1.253314*mad(x2,0));  % mean absolute deviation
end

end