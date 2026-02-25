clc; clear; close all;

% --- 参数设置 ---
randn('seed', 0);
rand('seed', 0);
K = 3;
M_x = 16;
M_y = 16;
T = 10;
c = 3e8;
f = 3e9;
lambda = c/f;
d = lambda/2;
Pos_x = d*(0:M_x-1).';
Pos_y = d*(0:M_y-1).';
theta = [30, 45, 60];
phi = [45, 75, 120];
SNR_list = 0:5:20; % 修改了一下范围起点，保证legend显示正常
numTrials = 200;

% --- 构造真实导向矢量 ---
A_x = zeros(M_x, K);
A_y = zeros(M_y, K);
for k = 1:K
    % 注意：这里定义了真实的物理模型
    % X轴对应 sin(theta)*sin(phi)
    A_x(:,k) = exp(1j*2*pi/lambda*Pos_x*sind(theta(k))*sind(phi(k)));
    % Y轴对应 sin(theta)*cos(phi)
    A_y(:,k) = exp(1j*2*pi/lambda*Pos_y*sind(theta(k))*cosd(phi(k)));
end

% 定义搜索网格 (传给搜索函数用)
search_grid_course = 0:1:180; 

h = waitbar(0, '开始蒙特卡洛仿真...');

for idx = 1:length(SNR_list)
    SNR = SNR_list(idx);
    mse_theta_sum = 0;
    mse_phi_sum   = 0;
    
   for trial=1:numTrials
        waitbar((idx-1)/length(SNR_list) + (trial/numTrials)*(1/length(SNR_list)), h, ...
            sprintf('SNR = %.1f dB | 试验 %d/%d | 总体进度: %.1f%%', ...
            SNR, trial, numTrials, 100*((idx-1)/length(SNR_list) + (trial/numTrials)*(1/length(SNR_list)))));
        % 信号生成
        S = randn(T, K) + 1j*randn(T, K);
        Y_clean = kr(A_y, A_x) * S.';
        Y = awgn(Y_clean, SNR, "measured");
        
        % 计算 CRB (保持原有逻辑不变)
        signal_power = mean(abs(Y_clean(:)).^2);
        noise_power = signal_power / (10^(SNR/10));
        X = S.';
        RS = (X*X')/T;
        AA = kr(A_y, A_x);
        R = (AA*X)*((AA*X)')/T + noise_power*eye(size(AA,1));
        
        % 计算导数矩阵用于 CRB
        A_x_thetad=zeros(size(A_x)); A_y_thetad=zeros(size(A_y));
        AA_theta=zeros(M_x*M_y,K);
        for k=1:K
            A_x_thetad(:,k)=(1j*2*pi/lambda*Pos_x*cosd(theta(k))*sind(phi(k))).*A_x(:,k);
            A_y_thetad(:,k)=(1j*2*pi/lambda*Pos_y*cosd(theta(k))*cosd(phi(k))).*A_y(:,k);
            AA_theta(:,k)=kron(A_y_thetad(:,k),A_x(:,k)) +kron(A_y(:,k),A_x_thetad(:,k));
        end
        
        A_x_phid=zeros(size(A_x)); A_y_phid=zeros(size(A_y));
        AA_phi=zeros(M_x*M_y,K);
        for k=1:K
            A_x_phid(:,k)=(1j*2*pi/lambda*Pos_x*sind(theta(k))*cosd(phi(k))).*A_x(:,k);
            A_y_phid(:,k)=(-1j*2*pi/lambda*Pos_y*sind(theta(k))*sind(phi(k))).*A_y(:,k);
            AA_phi(:,k)=kron(A_y_phid(:,k),A_x(:,k)) +kron(A_y(:,k),A_x_phid(:,k));
        end
        
        term1 = RS*AA'*inv(R);
        F_theta_theta=2*T*real((term1*AA_theta).*(term1*AA_theta).' + (term1*AA*RS).*(AA_theta'*inv(R)*AA_theta).');
        F_phi_phi=2*T*real((term1*AA_phi).*(term1*AA_phi).' + (term1*AA*RS).*(AA_phi'*inv(R)*AA_phi).');
        F_phi_theta=2*T*real((term1*AA_theta).*(term1*AA_phi).' + (term1*AA*RS).*(AA_theta'*inv(R)*AA_phi).');
        
        FIM=[F_theta_theta,F_phi_theta'; F_phi_theta,F_phi_phi];
        CRB_matrix=inv(FIM);
        crb_theta(trial,idx)=sum(diag(CRB_matrix(1:K, 1:K)))/K;
        crb_phi(trial,idx)=sum(diag(CRB_matrix(K+1:2*K, K+1:2*K)))/K;
        crb(trial,idx)=sum(diag(CRB_matrix))/(2*K);
        
        % --- 算法估计部分 (CPD + 1D Search) ---
        Y_tensor = reshape(Y, [M_y, M_x, T]);
        
        % 使用 CPD 分解 (tensorlab 或类似工具箱)
        % 注意：如果没有安装 tensorlab，这里会报错。
        % 假设已有 cpd 函数。
        U = cpd(Y_tensor, K); 
        A_y_hat = U{1}; % M_y * K
        A_x_hat = U{2}; % M_x * K
        
        est_theta = zeros(1, K);
        est_phi = zeros(1, K);
        
        for k = 1:K
            % 1. 对 X 轴导向矢量进行搜索
            % 搜索目标：找到 angle_x 使得 exp(j*k*x*cos(angle_x)) 匹配 A_x_hat
            % 物理意义：cos(angle_x) = sin(theta)*sin(phi)
            angle_x = fun_1D_search_single_deg(A_x_hat(:, k), search_grid_course, Pos_x, lambda);
            u_val = cosd(angle_x); % u = sin(theta)*sin(phi)
            
            % 2. 对 Y 轴导向矢量进行搜索
            % 搜索目标：找到 angle_y 使得 exp(j*k*y*cos(angle_y)) 匹配 A_y_hat
            % 物理意义：cos(angle_y) = sin(theta)*cos(phi)
            angle_y = fun_1D_search_single_deg(A_y_hat(:, k), search_grid_course, Pos_y, lambda);
            v_val = cosd(angle_y); % v = sin(theta)*cos(phi)
            
            % 3. 联合解算 Theta 和 Phi
            % 计算 Phi: tan(phi) = u/v -> phi = atan2(u, v)
            phi_est_k = atan2d(v_val, u_val);
            if phi_est_k < 0
                phi_est_k = phi_est_k + 360;
            end
            
            % 计算 Theta: u^2 + v^2 = sin^2(theta)*(sin^2 + cos^2) = sin^2(theta)
            sin_theta_sq = u_val^2 + v_val^2;
            % 限制数值误差防止 asin 出现复数
            if sin_theta_sq > 1
                sin_theta_sq = 1;
            end
            theta_est_k = asind(sqrt(sin_theta_sq));
            
            est_theta(k) = theta_est_k;
            est_phi(k)   = phi_est_k;
        end
        
        % 排序以匹配真实角度 (解决排列模糊)
        [est_theta_sorted, sort_index] = sort(est_theta);
        est_phi_sorted = est_phi(sort_index);
        
        mse_theta_sum = mse_theta_sum + sum((est_theta_sorted - theta).^2);
        mse_phi_sum   = mse_phi_sum   + sum((est_phi_sorted - phi).^2);
    end
    RMSE_theta(idx) = sqrt(mse_theta_sum / (numTrials * K));
    RMSE_phi(idx)   = sqrt(mse_phi_sum / (numTrials * K));
end

% --- 结果处理与绘图 ---
CRB_theta = sqrt((mean(crb_theta,1)));
CRB_phi   = sqrt((mean(crb_phi,1)));
CRB_total = sqrt((mean(crb,1)));
close(h);

figure;
% 注意：CRB计算出来是弧度，需要转角度显示
semilogy(SNR_list, CRB_theta*180/pi, '-o', 'LineWidth', 1.5); hold on;
semilogy(SNR_list, CRB_phi*180/pi,   '-*', 'LineWidth', 1.5); hold on;
semilogy(SNR_list, CRB_total*180/pi, '-^', 'LineWidth', 1.5); hold on;
semilogy(SNR_list, RMSE_theta, '--o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_list, RMSE_phi,   '--*', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
grid on;
xlabel('SNR (dB)');
ylabel('RMSE (Degrees)');
legend('CRB \theta', 'CRB \phi', 'CRB Total', 'RMSE \theta (Proposed)', 'RMSE \phi (Proposed)');
title('DOA Estimation Performance: Tensor Decomposition + 1D Search');

