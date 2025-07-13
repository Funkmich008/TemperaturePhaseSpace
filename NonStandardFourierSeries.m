function [a, b, g, AIC, BIC] = NonStandardFourierSeries(x, t, T, M, lambda)
% NonStandardFourierSeries
%   Fits a regularized truncated Fourier series of order M to the data (t, x)
%   with fundamental period T and regularization parameter lambda.
%   Returns:
%       a, b - Fourier cosine and sine coefficients
%       g    - function handle of the fitted Fourier series
%       AIC, BIC - information criteria for model selection
%
%   INPUTS:
%       x      - data vector of temperature anomalies (Nx1)
%       t      - time vector (Nx1)
%       T      - fundamental period of the Fourier series
%       M      - order of the Fourier series (number of harmonics)
%       lambda - regularization parameter controlling smoothness
%
%   OUTPUTS:
%       a, b   - vectors of cosine and sine coefficients
%       g      - fitted Fourier series as function handle g(t)
%       AIC, BIC - Akaike and Bayesian Information Criteria for fit quality

    % Ensure column vectors
    if size(x,2) > size(x,1), x = x.'; end
    if size(t,2) > size(t,1), t = t.'; end
    
    N = length(t);        % number of data points
    k = 0:M;              % harmonic indices
    L = length(k);        % number of cosine terms

    % Define angular frequencies for cosine and sine terms
    wc = 2*pi*(0:M)/T;    % frequencies for cosine terms
    ws = 2*pi*(1:M)/T;    % frequencies for sine terms

    % Compute right-hand side vector (projection of data onto Fourier basis)
    u_c = sum(x .* cos(wc.*t)).';    % projections onto cosines
    u_s = sum(x .* sin(ws.*t)).';    % projections onto sines
    u = [u_c; u_s];                  % system vector

    % Initialize design matrices
    Q = zeros(2*L-1, 2*L-1);   % data fidelity matrix
    R = zeros(2*L-1, 2*L-1);   % smoothness regularization matrix

    % Build Q and R matrices by summing over all data points
    for i = 1:N
        C = cos(wc*t(i));           % cosines at time t(i)
        S = sin(ws*t(i));           % sines at time t(i)
        dC2 = - wc.^2 .* C;         % second derivatives of cosines
        dS2 = - ws.^2 .* S;         % second derivatives of sines

        % Update data fidelity matrix Q
        Q = Q + [C.' * C, C.' * S; ...
                 S.' * C, S.' * S];

        % Update smoothness penalty matrix R
        R = R + [dC2.' * dC2, dC2.' * dS2; ...
                 dS2.' * dC2, dS2.' * dS2];
    end

    % Regularized normal equations matrix
    G = Q + lambda * R;

    % Solve for Fourier coefficients
    A = G \ u;                      % least-squares solution

    % Extract cosine (a) and sine (b) coefficients
    a = A(1:L);
    b = [0; A(L+1:2*L-1)];          % prepend b0=0 for indexing consistency

    % Construct function handle for fitted Fourier series
    g = @(xq) a(1) * ones(size(xq));
    for kidx = 2:L        
       w_k = 2*pi*(kidx-1)/T;
       g = @(xq) g(xq) + b(kidx) * sin(w_k * xq) + a(kidx) * cos(w_k * xq); 
    end

    % Compute residuals and information criteria
    residual = x - g(t);
    MSE = mean(residual.^2);
    K = length(A);                    % total number of fitted parameters
    AIC = N * log(MSE) + 2 * K;       % Akaike Information Criterion
    BIC = N * log(MSE) + K * log(N);  % Bayesian Information Criterion
    
end
