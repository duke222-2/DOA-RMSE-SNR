clc;clear;close all;
randn('seed', 0);
rand('seed', 0);
K=3;
M_x=16;
M_y=16;
T=100;
c=3e8;
f=3e9;
lambda=c/f;
d=lambda/2;
Pos_x=d*(0:M_x-1);
Pos_x=Pos_x.';
Pos_y=d*(0:M_y-1);
Pos_y=Pos_y.';
theta=[30,45,60];
phi=[45,75,120];
grid_theta=0:1:180;
grid_phi=0:1:180;
theta_tensor=zeros(1,K);
phi_tensor=zeros(1,K);
SNR_list=-0:5:20;
numTrials=200;


A_x=zeros(M_x,K);
A_y=zeros(M_y,K);
for k=1:K
    A_x(:,k)=exp(1j*2*pi/lambda*Pos_x*sind(theta(k))*sind(phi(k)));
    A_y(:,k)=exp(1j*2*pi/lambda*Pos_y*sind(theta(k))*cosd(phi(k)));
end
h = waitbar(0, '开始蒙特卡洛仿真...');
for idx=1:length(SNR_list)
    SNR=SNR_list(idx);
    mse_theta_sum = 0;
    mse_phi_sum   = 0;
    for trial=1:numTrials
        waitbar((idx-1)/length(SNR_list) + (trial/numTrials)*(1/length(SNR_list)), h, ...
            sprintf('SNR = %.1f dB | 试验 %d/%d | 总体进度: %.1f%%', ...
            SNR, trial, numTrials, 100*((idx-1)/length(SNR_list) + (trial/numTrials)*(1/length(SNR_list)))));
        S = randn(T,K) + 1j*randn(T,K);
        %{
S_power = zeros(1, K);
for k = 1:K
    S_power(k) = mean(abs(S(:, k)).^2);  % 计算第k个目标的平均功率
end
        %}
        Y_clean = kr(A_y, A_x) * S.';
        Y=awgn(Y_clean,SNR,"measured");
        signal_power = mean(abs(Y_clean(:)).^2);
        noise_power = signal_power / (10^(SNR/10));
        %noise_power_actual = mean(abs(Y(:) - Y_clean(:)).^2);
        Y_tensor=reshape(Y,[M_y,M_x,T]);
        U=cpd(Y_tensor,K);
        A_y_hat=U{1};
        A_x_hat=U{2};
        est_theta = zeros(1, K);
        est_phi = zeros(1, K);
        for k = 1:K
            % === 1. 相除提取相位 (这一步必须由 angle 算出弧度) ===
            % 这一步无法避免弧度，因为复数的相位就是弧度定义的

            % X轴处理
            vec_x = A_x_hat(:, k);
            ratio_x = mean(vec_x(2:end) ./ vec_x(1:end-1));
            psi_x_rad = angle(ratio_x); % 算出的是弧度 [-pi, pi]

            % Y轴处理
            vec_y = A_y_hat(:, k);
            ratio_y = mean(vec_y(2:end) ./ vec_y(1:end-1));
            psi_y_rad = angle(ratio_y); % 算出的是弧度 [-pi, pi]

            % === 2. 转换中间变量 (去量纲) ===
            % 你的阵元间距 d = lambda/2，导致相位系数是 pi
            % u, v 实际上代表的是 sin(theta)sin(phi) 和 sin(theta)cos(phi) 这种数值
            u = psi_x_rad / pi;
            v = psi_y_rad / pi;

            % === 3. 直接解算出"度数" ===

            % 使用 atan2d 直接得到度数
            % 注意顺序：atan2d(Y坐标也就是sin项, X坐标也就是cos项) -> (u, v)
            phi_degree = atan2d(v, u);

            % 修正负角度到 [0, 360]
            if phi_degree < 0
                phi_degree = phi_degree + 360;
            end

            % 计算 sin(theta) 的值
            sin_theta = sqrt(u^2 + v^2);

            % 保护一下防止超过 1
            if sin_theta > 1
                sin_theta = 1;
            end

            % 使用 asind 直接得到度数
            theta_degree = asind(sin_theta);

            est_theta(k) = theta_degree;
            est_phi(k)   = phi_degree;
        end
        [est_theta_sorted, sort_index] = sort(est_theta);
        est_phi_sorted = est_phi(sort_index);
        mse_theta_sum = mse_theta_sum + sum((est_theta_sorted - theta).^2);
        mse_phi_sum   = mse_phi_sum   + sum((est_phi_sorted - phi).^2);

    end
    RMSE_theta(idx) = sqrt(mse_theta_sum / (numTrials * K));
    RMSE_phi(idx)   = sqrt(mse_phi_sum / (numTrials * K));
end


close(h)
save('tensor2D.mat', 'RMSE_theta', 'RMSE_phi', 'SNR_list');
figure
semilogy(SNR_list, RMSE_theta, '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_list, RMSE_phi,   '-*', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;

grid on;
xlabel('SNR (dB)');
ylabel('RMSE (Degrees)');
legend('Theta RMSE', 'Phi RMSE');
title('2D DOA Estimation Performance');
ylim([10^(-3), 10^(2)]);
