addpath('C:\Users\Administrator\OneDrive\Semieff_CPL20250424\code')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Frank
N0 = 500; N = 50;
N1 = 10;
para = linspace(2,8,N1);
ple = zeros(N1,4); ple1 = zeros(N1,4);
kstep = zeros(N1,4); kstep1 = zeros(N1,4); kstep2 = zeros(N1,4); kstep3 = zeros(N1,4);
b1 = zeros(N0,1); b2 = zeros(N0,1); b3 = zeros(N0,1); b4 = zeros(N0,1);
LB = zeros(N1,1); g2 = zeros(N0,1);
for j = 1:N1
    alpha = para(j);
    z0 = 'F';
    LB(j) = Information_matrix(z0,alpha,1);
end

for j = 1:N1
    alpha = para(j);
    z0 = 'Frank';
    tic
parfor k = 1:N0
    rng(k,'twister')
    U1 = copularnd(z0,alpha,N);
    [~,I1] = sort(U1(:,1)); u = U1(:,1); v = U1(:,2);
    [~,I2] = sort(U1(:,2));
    [~,R1] = sort(I1);
    [~,R2] = sort(I2);
    R1 = R1/(N+1);
    R2 = R2/(N+1);
    U2 = [R1,R2];
    PML = copulafit(z0,U2);

    onestep = OSE(z0,1,PML,R1,R2,0,9,100); % 1 15 100 clayton Frank c - 0.5I(alpha>7)
    b1(k) = onestep;
    N3 = 9*(PML<5) + 12*(PML>5);
    b2(k) = OSE(z0,10,PML,R1,R2,0.5,N3,100); % j=6 c=0.5 N=12 0.87
    % b3(k) = OSE(z0,6,alpha0,R1,R2,0,N2,100);
    g2(k) = PML;
end
toc
std(b1)/std(g2) % 14 0.9852
std(b2)/std(g2)
kstep(j,1) = std(b1);
kstep(j,2) = mean((b1-alpha).*(b1-alpha));
kstep(j,3) = quantile(b1,0.7);
kstep(j,4) = mean(b1-alpha);
kstep2(j,1) = std(b2);
kstep2(j,2) = mean((b2-alpha).*(b2-alpha));
kstep2(j,3) = quantile(b2,0.7);
kstep2(j,4) = mean(b2-alpha);
kstep3(j,1) = std(b3);
kstep3(j,2) = mean((b3-alpha).*(b3-alpha));
kstep3(j,3) = quantile(b3,0.7);
kstep3(j,4) = mean(b3-alpha);
ple(j,1) = std(g2);
ple(j,2) = mean((g2-alpha).*(g2-alpha));
ple(j,3) = quantile(g2,0.5);
ple(j,4) = mean(g2-alpha);
ple1(j,1) = std(p1);
ple1(j,2) = mean((p1-alpha).*(p1-alpha));
ple1(j,3) = quantile(p1,0.5);
ple1(j,4) = mean(p1-alpha);
end

save('kstep_frank.mat', 'kstep'); 
save('kstep_frank2.mat', 'kstep2'); 
save('plm_frank.mat', 'ple'); 
% Plotting
% kstep = load('kstep_frank.mat');
% kstep = kstep.kstep;
% ple = load('plm_frank.mat');
% ple = ple.ple;
v = 1;
figure;
hold on;
plot(para, kstep(:,v), 'b-o', 'LineWidth', 2, 'DisplayName', 'K=1');
plot(para, kstep2(:,v), 'c-s', 'LineWidth', 2, 'DisplayName', 'K=2');
plot(para, ple(:,v), 'r--s', 'LineWidth', 2, 'DisplayName', 'PLE');
plot(para,power(N*LB,-1/2), 'b', 'LineWidth', 2, 'DisplayName', 'OSE');

% Use LaTeX format for θ symbol
xlabel('\theta', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Standard Variance', 'FontSize', 12);
title('Standard Variance Comparison', 'FontSize', 14);
legend('Location', 'northwest');
grid on;
hold off;


v = 2;
figure;
hold on;
plot(para, kstep(:,v), 'b-o', 'LineWidth', 2, 'DisplayName', 'kstep');
plot(para, ple(:,v), 'r--s', 'LineWidth', 2, 'DisplayName', 'ple');
% plot(para, power(N*LB,-1/2), 'g:', 'LineWidth', 2, 'DisplayName', 'LB');

% Place legend in the upper-left corner
legend('Location', 'northwest'); % Key modification

xlabel('\theta');
ylabel('MSE');
title('MSE Comparison');
grid on;
hold off;


% save('kstep_frank.mat', 'kstep'); 
% save('plm_frank.mat', 'ple'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Clayton 
N0 = 500; N = 50;
N1 = 10;
para = linspace(2,10,N1);
ple = zeros(N1,4); ple1 = zeros(N1,4);
kstep = zeros(N1,4); kstep1 = zeros(N1,4); kstep2 = zeros(N1,4); kstep3 = zeros(N1,4);
b1 = zeros(N0,1); b2 = zeros(N0,1); b3 = zeros(N0,1); b4 = zeros(N0,1);
eff_vec = zeros(N1,1);
g2 = zeros(N0,1);
for j = 1:N1
    alpha = para(j);
    z0 = 'C';
    LB(j) = Information_matrix(z0,alpha,1);
end

for j = 1:N1
    alpha = para(j);
    tic
parfor k = 1:N0
    rng(k, 'twister')
    z0 = 'Clayton';
    U1 = copularnd(z0,alpha,N);
    [~,I1] = sort(U1(:,1)); u = U1(:,1); v = U1(:,2);
    [~,I2] = sort(U1(:,2));
    [~,R1] = sort(I1);
    [~,R2] = sort(I2);
    R1 = R1/(N+1);
    R2 = R2/(N+1);
    U2 = [R1,R2];
    PML = copulafit(z0,U2);
    B = 100;
    AB_boot = zeros(B,1); AB_boot1 = zeros(B,1);
    for jj = 1:B
        boot_idx = randsample(N, N, true);
        U2_boot = U2(boot_idx, :);
        PML = copulafit(z0,U2_boot);
        AB_boot(jj) = OSE(z0, 1, PML, U2_boot(:,1), U2_boot(:,2), 0, 14, 100);
        AB_boot1(jj) = OSE(z0, 1, PML, U2_boot(:,1), U2_boot(:,2), 0,9, 100);
    end
    if std(AB_boot) > std(AB_boot1)
        N2 = 9;
    else
        N2 = 14;
    end
    PML;
    AB = OSE(z0,1,PML,R1,R2,0,N2,100);
    b1(k) = AB;
    g2(k) = PML;
end
toc
std(b1)/std(g2) % 14 0.9852
std(b2)/std(g2)
std(b3)/std(g2)
kstep(j,1) = std(b1);
kstep(j,2) = mean((b1-alpha).*(b1-alpha));
kstep(j,3) = quantile(b1,0.7);
kstep(j,4) = mean(b1-alpha);
kstep2(j,1) = std(b2);
kstep2(j,2) = mean((b2-alpha).*(b2-alpha));
kstep2(j,3) = quantile(b2,0.7);
kstep2(j,4) = mean(b2-alpha);
kstep3(j,1) = std(b3);
kstep3(j,2) = mean((b3-alpha).*(b3-alpha));
kstep3(j,3) = quantile(b3,0.7);
kstep3(j,4) = mean(b3-alpha);
ple(j,1) = std(g2);
ple(j,2) = mean((g2-alpha).*(g2-alpha));
ple(j,3) = quantile(g2,0.5);
ple(j,4) = mean(g2-alpha);
ple1(j,1) = std(p1);
ple1(j,2) = mean((p1-alpha).*(p1-alpha));
ple1(j,3) = quantile(p1,0.5);
ple1(j,4) = mean(p1-alpha);
end
find(kstep3(:,1)-ple(:,1)>0)

% save('kstep_clayton.mat', 'kstep'); 
% save('plm_clayton.mat', 'ple'); 
% Plotting
% kstep = load('kstep_clayton.mat');
% kstep = kstep.kstep;
% ple = load('plm_clayton.mat');
% ple = ple.ple;
v = 2;
figure;
hold on;
plot(para, kstep(:,v), 'b-o', 'LineWidth', 2, 'DisplayName', 'OSE');
plot(para, ple(:,v), 'r--s', 'LineWidth', 2, 'DisplayName', 'PLE');


xlabel('\theta');
ylabel('MSE', 'FontSize', 12);
title('MSE Comparison', 'FontSize', 14);
legend('Location', 'northwest');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% % t
N0 = 500; N = 100;
N1 = 10;
N1 = 10; kstep = zeros(N1,4); kstep2 = zeros(N1,4);
ple = zeros(N1,4); ose1 = zeros(N1,4);
b1 = zeros(N0,1); b2 = zeros(N0,1); b3 = zeros(N0,1);
para = linspace(-0.9,0.9,N1);
nu = 5; N = 50;
for j = 1:N1
    alpha = para(j);
    z0 = 't';
    LB(j) = Information_matrix(z0,alpha,nu);
end

for j = 9:N1
    alpha = para(j);
    z0 = 't';
parfor k = 1:N0
    U1 = copularnd(z0,alpha,nu,N);
    [~,I1] = sort(U1(:,1)); u = U1(:,1); v = U1(:,2);
    [~,I2] = sort(U1(:,2));
    [~,R1] = sort(I1);
    [~,R2] = sort(I2);
    R1 = R1/(N+1);
    R2 = R2/(N+1);
    U2 = [R1,R2];
    % Assume U2 is an N×2 matrix
    NZ = size(U2, 1); % Get total sample size
    k0 = 40; % Number of samples to draw
    % Randomly draw 40 unique indices
    idx = randperm(NZ, k0); 
    % Extract corresponding samples
    U2_sampled = U2(idx, :);
    tau = corr(U2_sampled(:,1), U2_sampled(:,2), 'Type', 'Kendall');
    Rho_initial = sin(pi * tau / 2);
    Rho_initial(Rho_initial >= 1) = 0.99;
    Rho_initial(Rho_initial <= -1) = -0.99;
    onestep = OSEt(1,Rho_initial,5,R1,R2,20);
    b1(k) = onestep;
    b2(k) = Rho_initial;
    b3(k) = OSEt(10,Rho_initial,5,R1,R2,20);
end
kstep(j,1) = std(b1);
kstep(j,2) = mean((b1-alpha).*(b1-alpha));
kstep(j,3) = quantile(b1,0.5);
kstep(j,4) = mean(b1-alpha);
kstep2(j,1) = std(b3);
kstep2(j,2) = mean((b3-alpha).*(b3-alpha));
kstep2(j,3) = quantile(b3,0.5);
kstep2(j,4) = mean(b3-alpha);
ple(j,1) = std(b2);
ple(j,2) = mean((b2-alpha).*(b2-alpha));
ple(j,3) = quantile(b2,0.5);
ple(j,4) = mean(b2-alpha);
end

% save('kstep_t.mat', 'kstep'); 
% save('kstep2_t.mat', 'kstep2'); 
% save('plm_t.mat', 'ple'); 
save('kstep_t_n100.mat', 'kstep'); 
save('kstep2_t_n100.mat', 'kstep2'); 
save('plm_t.mat_n100', 'ple'); 
% Plotting
kstep = load('kstep_t.mat');
kstep = kstep.kstep;
kstep2 = load('kstep2_t.mat');
kstep2 = kstep2.kstep2;
ple = load('plm_t.mat');
ple = ple.ple;
v = 1;
% Clear figure window and set global properties
figure;
hold on;

% Use different line styles, markers, and colors
p2 = plot(para, ple(:,v), 'r--s', 'LineWidth', 2, 'MarkerSize', 6, ...
          'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'K=0');
      
p1 = plot(para, kstep(:,v), 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
          'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'K=1');
      
p3 = plot(para, kstep2(:,v), 'm-.^', 'LineWidth', 2, 'MarkerSize', 6, ...
          'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'DisplayName', 'K=2');
      
p4 = plot(para, power(N*LB,-1/2), 'k:d', 'LineWidth', 2, 'MarkerSize', 6, ...
          'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'DisplayName', 'SEB');

% Fix legend position setting (MATLAB 2016 compatible)
legend('K=0','K=1', 'K=2','SEB', 10, 'Box', 'off');

% Set axis labels
xlabel('\theta', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Standard Variance', 'FontSize', 12, 'FontWeight', 'bold');

% Grid settings
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

% Set figure background color
set(gcf, 'Color', [0.95 0.95 0.95]);

% Adjust axis limits
xlim([min(para) max(para)]);
ylim auto;

% Add title
title('Comparison of Different Methods', 'FontSize', 14, 'FontWeight', 'bold');

% Fix transparency setting (MATLAB 2016 compatible)
set(p1, 'Color', [0 0 1 0.7]); % Blue with transparency
set(p2, 'Color', [1 0 0 0.7]); % Red with transparency
set(p3, 'Color', [1 0 1 0.7]); % Magenta with transparency
set(p4, 'Color', [0 0 0 0.7]); % Black with transparency