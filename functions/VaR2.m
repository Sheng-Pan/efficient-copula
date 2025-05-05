function [q] = VaR2(copula_type, returns, theta, confidence_level)
% Parameter settings
M = 100000; % Number of simulated samples
w = [0.5, 0.5]; % Portfolio weights

% Simulated data (or replace with actual data)

% Clean data (remove NaN and Inf)
returns = returns(all(~isnan(returns) & ~isinf(returns), 2), :);
N = size(returns, 1);

% Step 1: Transform to uniform distribution
U = zeros(N, 2);
for i = 1:2
    [F, X] = ecdf(returns(:,i));
    [X, idx] = unique(X, 'stable'); % Remove duplicates
    F = F(idx);
    U(:,i) = interp1(X, F, returns(:,i), 'linear', 'extrap');
    U(:,i) = max(min(U(:,i), 1-eps), eps); % Restrict to (0,1)
end

% Step 3: Generate simulated data
U_sim = copularnd(copula_type, theta, M);
U_sim = max(min(U_sim, 1-eps), eps); % Restrict U_sim
returns_sim = zeros(M, 2);
for i = 1:2
    [F, X] = ecdf(returns(:,i));
    [X, idx] = unique(X, 'stable');
    F = F(idx);
    returns_sim(:,i) = interp1(F, X, U_sim(:,i), 'linear', 'extrap');
end

% Step 4: Calculate portfolio returns and losses
portfolio_returns = w(1) * returns_sim(:,1) + w(2) * returns_sim(:,2);
losses = -portfolio_returns;

% Step 5: Estimate VaR
q = quantile(losses, 1 - confidence_level);
end