function [Xgrid, Ygrid, pdf] = ComputePhaseSpaceDensity(X, Y, K)
% ComputePhaseSpaceDensity
%   Estimates the probability density function (PDF) in the 2D phase space
%   defined by X and Y, using kernel density estimation (ksdensity).
%   Returns:
%       Xgrid, Ygrid - evaluation grids for plotting
%       pdf          - estimated density on the grid
%
%   INPUTS:
%       X, Y - vectors or matrices of phase space coordinates
%       K    - number of grid points per axis (resolution)
%
%   OUTPUTS:
%       Xgrid, Ygrid - 2D grids of coordinates (KxK)
%       pdf          - 2D matrix of estimated density values (KxK)

    % Combine X and Y into a single point cloud for 2D kernel density estimation
    xy = [X(:), Y(:)];
    
    % Define expanded support for density estimation (15% padding on each side)
    xmin = min(X(:)) - 0.15 * abs(min(X(:)));
    xmax = max(X(:)) + 0.15 * abs(max(X(:)));
    ymin = min(Y(:)) - 0.15 * abs(min(Y(:)));
    ymax = max(Y(:)) + 0.15 * abs(max(Y(:)));
    support = [xmin ymin; xmax ymax];

    % Generate uniform grid on the extended support
    Tx = linspace(xmin, xmax, K);
    Vy = linspace(ymin, ymax, K);
    [Xgrid, Ygrid] = meshgrid(Tx, Vy);
    
    % Estimate density at each grid point
    pdf = ksdensity(xy, [Xgrid(:), Ygrid(:)], 'Support', support);
    pdf = reshape(pdf, size(Xgrid));
end
