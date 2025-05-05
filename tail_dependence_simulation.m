N1 = 10;
alpha_values = linspace(2,10,N1);
N0 = 500;            % Number of repeated experiments
N = 300;            % Sample size
B = 10;             % 10-fold cross-validation
q = 0.95;           % Tail dependence quantile threshold

z0 = 'Frank';       % Copula type used for data generation
z1 = 'Clayton';     % Copula type used for fitting
results = struct();

for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    N2 = 14*(a_idx<=4) + 9*(a_idx>4);

    % Pre-allocate storage variables
    TD_vec = zeros(N0, 1);      % Store final selected tail dependence estimates
    TD_mis = zeros(N0, 1);      % Store alternative tail dependence estimates
    mse_nonpa_all = zeros(N0, 1); % MSE for non-parametric method
    mse_ose_all = zeros(N0, 1);   % MSE for OSE method
    wk = zeros(N0, 1);
    parfor k = 1:N0  % Parallelize to speed up repeated experiments
        rng(k, 'twister');  % Fix random seed
        % Generate Copula data (Frank)
        U = copularnd(z0, alpha, N);
        
        % Create cross-validation partitions
        cv = cvpartition(N, 'KFold', B);
        
        % Initialize MSE storage
        mse_nonpa = zeros(B, 1);
        mse_ose = zeros(B, 1);
        
        for fold = 1:B
            % Split data
            trainIdx = cv.training(fold);
            testIdx = cv.test(fold);
            U_train = U(trainIdx, :);
            U_test = U(testIdx, :);
              
            R1_train = tiedrank(U_train(:,1)) / (length(U_train(:,1)) + 1);
            R2_train = tiedrank(U_train(:,2)) / (length(U_train(:,1)) + 1);
              
            R1_test = tiedrank(U_test(:,1)) / (length(U_test(:,1)) + 1);
            R2_test = tiedrank(U_test(:,2)) / (length(U_test(:,1)) + 1);
            % --- Non-parametric tail dependence estimation ---
            % Lower tail dependence
            lower_tail = (R1_train <= (1-q)) & (R2_train <= (1-q));
            TD_nonpa = sum(lower_tail) / (length(U_train(:,1)) * (1 - q));
            % --- OSE method tail dependence estimation ---
 
            % Fit Clayton Copula
            param = copulafit(z1, U2_train);
            param_ose = OSE(z1, 1, param, R1_train, R2_train, 0,N2,100);
            
            % Monte Carlo simulation for tail dependence
            TD_ose = 2^(-1/param_ose);
            
            % True tail dependence on test set
            lower_tail = (R1_test <= (1-q)) & (R2_test <= (1-q));
            TD_test = sum(lower_tail) / (length(U_test(:,1)) * (1 - q));
            
            % Record MSE
            mse_nonpa(fold) = (TD_test - TD_nonpa)^2;
            mse_ose(fold) = (TD_test - TD_ose)^2;
        end
        
        % --- Tail dependence estimation on full dataset ---
        % Method 1: Non-parametric
        R1 = tiedrank(U(:,1)) / (length(U(:,1)) + 1);
        R2 = tiedrank(U(:,2)) / (length(U(:,1)) + 1);
        lower_tail = (R1 <= (1-q)) & (R2 <= (1-q));
        TD_nonpa_full = sum(lower_tail) / (length(U(:,1)) * (1 - q));
        
        % Method 2: OSE method
     
        param_full = copulafit(z1, [R1, R2]);
        param_ose = OSE(z1, 1, param, R1, R2, 0,N2,100);
        TD_ose_full = 2^(-1/param_ose);
        
        % Select estimation method
        if mean(mse_nonpa) >= mean(mse_ose)
            TD_vec(k) = TD_ose_full;
        else
            TD_vec(k) = TD_nonpa_full;
        end
        wk(k) = mean(mse_nonpa) >= mean(mse_ose);
        % PLE
        TD_mis(k) = 2^(-1/param_full); % sum(lower_tail) / (length(U_sim(:,1)) * (1 - q));
    end
    
    % --- Calculate true tail dependence (theoretical value) ---
    % Theoretical upper tail dependence for Clayton Copula
    TD_true = 0; % 2^(-1/alpha); % Formula for Clayton's lower tail dependence
    
    % Store results
    results(a_idx).alpha = alpha;
    results(a_idx).mse_vec = mean((TD_vec - TD_true).^2);
    results(a_idx).mse_mis = mean((TD_mis - TD_true).^2);
    results(a_idx).mean_TD_vec = mean(TD_vec);
    results(a_idx).mean_TD_mis = mean(TD_mis);
    results(a_idx).TD_true = TD_true;
    fprintf('\r[', a_idx, '/', N, '] '); % Return to start of line
end
% save('Frank_Clayton.mat', 'results'); 
rr = load('Frank_Clayton.mat');
results = rr.results;
N1 = 10;
alpha_values = linspace(2,10,N1);

% Display results
disp('Alpha | MSE(Vec) | MSE(Mis) | Mean Vec | Mean Mis | True TD');
for a_idx = 1:length(alpha_values)
    fprintf('%5.2f | %9.4f | %9.4f | %8.4f | %8.4f | %8.4f\n', ...
            results(a_idx).alpha, ...
            results(a_idx).mse_vec, ...
            results(a_idx).mse_mis, ...
            results(a_idx).mean_TD_vec, ...
            results(a_idx).mean_TD_mis, ...
            results(a_idx).TD_true);
end

% Plotting
figure;
plot([results.alpha], [results.mse_vec], 'b-o', 'DisplayName', 'MV');
hold on;
plot([results.alpha], [results.mse_mis], 'r-s', 'DisplayName', 'PLE');
xlabel('\theta');
ylabel('MSE');
title('Tail Dependence Estimation: Misspecification of Frank Copula as Clayton', ...
    'FontSize', 12); % Set title font size to 12
legend('Location', 'northwest');
grid on;