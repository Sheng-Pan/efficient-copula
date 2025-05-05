N1 = 10;
alpha_values = linspace(2,10,N1);
N0 = 500;            % Number of repeated experiments
N = 200;            % Sample size
B = 5;              % 5-fold cross-validation
conf_level = 0.95;  % VaR confidence level

z0 = 'Gumbel';     % Copula type used for data generation
z1 = 'Clayton';    % Copula type used for fitting
results = struct();

for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    N2 = 14*(a_idx<=4) + 9*(a_idx>4);

    % Pre-allocate storage variables
    VaR_vec = zeros(N0, 1);      % Store final selected VaR estimates
    VaR_mis = zeros(N0, 1);      % Store alternative VaR estimates
    mse_nonpa_all = zeros(N0, 1); % MSE for non-parametric method
    mse_ose_all = zeros(N0, 1);   % MSE for OSE method

    parfor k = 1:N0  % Parallelize to speed up repeated experiments
        rng(k, 'twister');  % Fix random seed
        % Generate Copula data (Gumbel)
        U1 = copularnd(z0, alpha, N);
        
        % Transform to asset returns (standard normal marginals)
        returns = norminv(U1, 0, 1);
        
        % Create cross-validation partitions
        cv = cvpartition(N, 'KFold', B);
        
        % Initialize MSE storage
        mse_nonpa = zeros(B, 1);
        mse_ose = zeros(B, 1);
        
        for fold = 1:B
            % Split data
            trainIdx = cv.training(fold);
            testIdx = cv.test(fold);
            trainReturns = returns(trainIdx, :);
            testReturns = returns(testIdx, :);
            
            % --- Non-parametric VaR calculation ---
            portfolio_train = 0.5*trainReturns(:,1) + 0.5*trainReturns(:,2);
            VaR_nonpa = quantile(-portfolio_train, 1-conf_level);
            
            % --- OSE method VaR calculation ---
            % Compute empirical CDF
            % Method 1: Non-parametric
            R1 = tiedrank(trainReturns(:,1)) / (length(trainReturns(:,1)) + 1);
            R2 = tiedrank(trainReturns(:,2)) / (length(trainReturns(:,1)) + 1);
            U2_train = [R1, R2];
            
            % Fit Clayton Copula
            param = copulafit(z1, U2_train);
            param_ose = max(OSE(z1, 1, param, R1, R2, 0,N2,100),1e-05);
            
            % Monte Carlo simulation for VaR
            VaR_ose = VaR2(z1,trainReturns,param_ose,0.95);
            
            % True VaR on test set
            portfolio_test = 0.5*testReturns(:,1) + 0.5*testReturns(:,2);
            VaR_test = quantile(-portfolio_test,1-conf_level);
            
            % Record MSE
            mse_nonpa(fold) = (VaR_test - VaR_nonpa)^2;
            mse_ose(fold) = (VaR_test - VaR_ose)^2;
        end
        
        % --- VaR estimation on full dataset ---
        % Method 1: Non-parametric
        full_portfolio = 0.5*returns(:,1) + 0.5*returns(:,2);
        VaR_nonpa_full = quantile(-full_portfolio, 1-conf_level);
        
        % Method 2: OSE method
        R1 = tiedrank(returns(:,1)) / (length(returns(:,1)) + 1);
        R2 = tiedrank(returns(:,2)) / (length(returns(:,1)) + 1);
        param_full = max(copulafit(z1, [R1, R2]),1e-6);
        param_ose = max(OSE(z1, 1, param_full, R1, R2, 0,N2,100),1e-6);
        VaR_ose_full = VaR2(z1,returns,param_ose,0.95);
        
        % Select VaR calculation method
        if mean(mse_nonpa) >= mean(mse_ose)
            VaR_vec(k) = VaR_ose_full;
        else
            VaR_vec(k) = VaR_nonpa_full;
        end
        
        % Alternative VaR estimate (directly using fitted Copula parameters)
        VaR_mis(k) = VaR2(z1,returns,param_full,0.95);
        
    end
    
    % --- Calculate true VaR (theoretical value) ---
    U_true = copularnd(z0, alpha, 10000);
    true_returns = norminv(U_true, 0, 1);
    VaR_true = quantile(-(0.5*true_returns(:,1)+0.5*true_returns(:,2)), 1-conf_level);

    % Store results
    results(a_idx).alpha = alpha;
    results(a_idx).mse_vec = mean((VaR_vec - VaR_true).^2);
    results(a_idx).mse_mis = mean((VaR_mis - VaR_true).^2);
    results(a_idx).mean_VaR_vec = mean(VaR_vec);
    results(a_idx).mean_VaR_mis = mean(VaR_mis);
    results(a_idx).VaR_true = VaR_true;
end

% save('Gumbel_Clayton.mat', 'results'); 
rr = load('Gumbel_Clayton.mat');
results = rr.results;
N1 = 10;
alpha_values = linspace(2,10,N1);

% Display results
disp('Alpha | MSE(Vec) | MSE(Mis) | Mean Vec | Mean Mis | True VaR');
for a_idx = 1:length(alpha_values)
    fprintf('%5.2f | %9.4f | %9.4f | %8.4f | %8.4f | %8.4f\n', ...
            results(a_idx).alpha, ...
            results(a_idx).mse_vec, ...
            results(a_idx).mse_mis, ...
            results(a_idx).mean_VaR_vec, ...
            results(a_idx).mean_VaR_mis, ...
            results(a_idx).VaR_true);
end

% Plotting
figure;
plot([results.alpha], [results.mse_vec], 'b-o', 'DisplayName', 'MV');
hold on;
plot([results.alpha], [results.mse_mis], 'r-s', 'DisplayName', 'PLE');
xlabel('\theta');
ylabel('MSE');
title('VaR Estimation: Misspecification of Gumbel Copula as Clayton', ...
    'FontSize', 12);
legend('Location', 'northwest');
grid on;