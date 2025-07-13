function [X_space, Xi_space, pdf_matrix] = ComputeDensity(X, Y, xiLim, K)
% ComputeDensity
%   Estimates time-resolved 1D probability densities (PDFs) for ensemble
%   samples Y at each time point X, using kernel density estimation (KDE).
%   Returns grids suitable for density carpet plots.
%
%   INPUTS:
%       X      - vector of time points (1 x nTime)
%       Y      - matrix of ensemble samples (nEnsemble x nTime)
%       xiLim  - two-element vector specifying min and max of xi grid [xiMin xiMax]
%       K      - number of points in xi grid
%
%   OUTPUTS:
%       X_space, Xi_space - grids for plotting (K x nTime)
%       pdf_matrix        - density matrix (K x nTime), where each column is the
%                           estimated PDF for the corresponding time point

    nTime = length(X);                  % number of time points
    xi = linspace(xiLim(1), xiLim(2), K);  % evaluation grid for density
    pdf_matrix = zeros(length(xi), nTime); % initialize density matrix
    
    for i = 1:nTime
        samples = Y(:, i);              % ensemble samples at time X(i)
        [f, ~] = ksdensity(samples, xi);% kernel density estimate at grid points xi
        pdf_matrix(:, i) = f;           % store PDF for current time point
        
        % Optional: display progress
        disp(['Progress: ', num2str(round(i/nTime*100, 2)), '%'])
    end

    % Create meshgrids for easy contour or surface plotting
    [X_space, Xi_space] = meshgrid(X, xi);
end
