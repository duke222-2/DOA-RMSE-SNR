clc;
clear;
tic

%% 入射信号参数设置
c = 3e8;
f = 3e9;
lambda = c/f;
theta =[20,75,150]; % 方位角（度）
phi = [20,75,150];   % 俯仰角（度）
K = length(theta);
x_True = sind(phi).*sind(theta);
y_True = cosd(phi).*sind(theta);

%% 平面阵列参数设置
N_x = 8;
N_y = 8;
d_x = lambda/2;
d_y = lambda/2;
L = 20; % 固定快拍数
SNR_list = 0:5:20; % SNR范围从-10dB到20dB，步长5dB

%% 互耦系数设置
kappa = 1;  %沿x轴方向的互耦系数的个数
gamma = 1;  %沿y轴方向的互耦系数的个数
c1 = [1; 0.3527+0.4854i; zeros(N_x-kappa-1,1)];
c2 = [0.3527+0.4854i; 0.0927-0.2853i; zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block = cell(N_y,1); %构造互耦数胞
C_block{1} = toeplitz(c1,c1);
C_block{2} = toeplitz(c2,c2);
for i = gamma + 2 : N_y
    C_block{i} = zeros(N_x,N_x);
end
%构造分块互耦系数矩阵
C = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block{idx};
    end
end

%% 接收矩阵的构造
A = zeros(N_x*N_y,K);
for k = 1:K
    a_x = exp(-1j*2*pi/lambda*d_x * (0:N_x-1)'.*x_True(k));
    a_y = exp(-1j*2*pi/lambda*d_y * (0:N_y-1)'.*y_True(k));
    A(:,k) = kron(a_y,a_x);
end
A_bar = C*A;

%% 构造A_dot_bar_theta
A_dot_theta = zeros(N_x*N_y,K);
for k = 1:K
    for ny = 1:N_y
        for nx = 1:N_x
            A_dot_theta(nx+(ny-1)*N_x,k) = ((-1j*2*pi/lambda*d_x*sind(phi(k))*cosd(theta(k ))*(nx-1))+(-1j*2*pi/lambda*d_y*cosd(phi(k))*cosd(theta(k ))*(ny-1)))*A(nx+(ny-1)*N_x,k);
        end
    end
end
A_dot_bar_theta = C*A_dot_theta;

%% 构造A_dot_bar_phi
A_dot_phi = zeros(N_x*N_y,K);
for k = 1:K
    for ny = 1:N_y
        for nx = 1:N_x
            A_dot_phi(nx+(ny-1)*N_x,k) = ((-1j*2*pi/lambda*d_x*cosd(phi(k))*sind(theta(k ))*(nx-1))+(1j*2*pi/lambda*d_y*sind(phi(k))*sind(theta(k ))*(ny-1)))*A(nx+(ny-1)*N_x,k);
        end
    end
end
A_dot_bar_phi = C*A_dot_phi;

%% 为对互耦矩阵求导作准备工作
% 对a1，b1求导的准备
c1_1a = [0; 1; zeros(N_x-kappa-1,1)];
c2_1a = [0; 0; zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block_1a = cell(N_y,1); %构造互耦数胞
C_block_1a{1} = toeplitz(c1_1a,c1_1a);
C_block_1a{2} = toeplitz(c2_1a,c2_1a);
for i = gamma + 2 : N_y
    C_block_1a{i} = zeros(N_x,N_x);
end
C_1a = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C_1a((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block_1a{idx};
    end
end

c1_1b = [0; sqrt(-1); zeros(N_x-kappa-1,1)];
c2_1b = [0; 0; zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block_1b = cell(N_y,1); %构造互耦数胞
C_block_1b{1} = toeplitz(c1_1b,c1_1b);
C_block_1b{2} = toeplitz(c2_1b,c2_1b);
for i = gamma + 2 : N_y
    C_block_1b{i} = zeros(N_x,N_x);
end
C_1b = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C_1b((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block_1b{idx};
    end
end

% 对a2，b2求导的准备
c1_2a = [0; 0; zeros(N_x-kappa-1,1)];
c2_2a = [1; 0; zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block_2a = cell(N_y,1); %构造互耦数胞
C_block_2a{1} = toeplitz(c1_2a,c1_2a);
C_block_2a{2} = toeplitz(c2_2a,c2_2a);
for i = kappa + 2 : N_y
    C_block_2a{i} = zeros(N_x,N_x);
end
C_2a = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C_2a((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block_2a{idx};
    end
end

c1_2b = [0; 0; zeros(N_x-kappa-1,1)];
c2_2b = [sqrt(-1); 0; zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block_2b = cell(N_y,1); %构造互耦数胞
C_block_2b{1} = toeplitz(c1_2b,c1_2b);
C_block_2b{2} = toeplitz(c2_2b,c2_2b);
for i = kappa + 2 : N_y
    C_block_2b{i} = zeros(N_x,N_x);
end
C_2b = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C_2b((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block_2b{idx};
    end
end

% 对a3，b3求导的准备
c1_3a = [0; 0; zeros(N_x-kappa-1,1)];
c2_3a = [0; 1; zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block_3a = cell(N_y,1); %构造互耦数胞
C_block_3a{1} = toeplitz(c1_3a,c1_3a);
C_block_3a{2} = toeplitz(c2_3a,c2_3a);
for i = kappa + 2 : N_y
    C_block_3a{i} = zeros(N_x,N_x);
end
C_3a = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C_3a((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block_3a{idx};
    end
end

c1_3b = [0; 0; zeros(N_x-kappa-1,1)];
c2_3b = [0; sqrt(-1); zeros(N_x-kappa-1,1)]; %非零互耦矩阵参数设置
C_block_3b = cell(N_y,1); %构造互耦数胞
C_block_3b{1} = toeplitz(c1_3b,c1_3b);
C_block_3b{2} = toeplitz(c2_3b,c2_3b);
for i = kappa + 2 : N_y
    C_block_3b{i} = zeros(N_x,N_x);
end
C_3b = zeros(N_x*N_y,N_x*N_y);
for i = 1:N_y
    for j = 1:N_y
        idx = abs(i-j)+1;
        C_3b((i-1)*N_x+1:i*N_x , (j-1)*N_x+1:j*N_x) = C_block_3b{idx};
    end
end

%% 构造A_dot_bar_a1,A_dot_bar_b1,A_dot_bar_a2,A_dot_bar_b2
A_dot_bar_a1 = C_1a*A;
A_dot_bar_b1 = C_1b*A;
A_dot_bar_a2 = C_2a*A;
A_dot_bar_b2 = C_2b*A;
A_dot_bar_a3 = C_3a*A;
A_dot_bar_b3 = C_3b*A;

%% 初始化CRB结果
CRB_DOA = zeros(length(SNR_list),1);
CRB_C = zeros(length(SNR_list),1);
norm_C = norm([1,0.3527+0.4854i,0.3527+0.4854i,0.0927-0.2853i],2);

%% 主循环：遍历不同的SNR
for snr_index = 1:length(SNR_list)
    SNR = SNR_list(snr_index);
    fprintf('计算SNR = %d dB...\n', SNR);
    
    % 将SNR(dB)转换为线性功率比
    Power = 10^(SNR/10);
    
    % 噪声功率设为1，信号功率根据SNR计算
    sigma_n = 1;  % 噪声方差设为1
    sigma_s = sigma_n * Power;  % 信号方差
    
    %% 计算协方差矩阵
    R_s = sigma_s * eye(K);    % 信号协方差矩阵
    R_x = C * A * R_s * A' * C' + sigma_n * eye(N_x*N_y);
    
    %% 构造F_theta_theta
    F_theta_theta = 2*L*real((R_s*(A_bar')*inv(R_x)*A_dot_bar_theta).*conj(R_s*(A_bar')*inv(R_x)*A_dot_bar_theta)' + (R_s*A_bar'*inv(R_x)*A_bar*R_s).*conj(A_dot_bar_theta'*inv(R_x)*A_dot_bar_theta)');
    
    %% 构造F_phi_phi
    F_phi_phi = 2*L*real((R_s*(A_bar')*inv(R_x)*A_dot_bar_phi).*conj(R_s*(A_bar')*inv(R_x)*A_dot_bar_phi)' + (R_s*A_bar'*inv(R_x)*A_bar*R_s).*conj(A_dot_bar_phi'*inv(R_x)*A_dot_bar_phi)');
    
    %% 构造F_phi_theta
    F_phi_theta = 2*L*real((R_s*(A_bar')*inv(R_x)*A_dot_bar_theta).*conj(R_s*(A_bar')*inv(R_x)*A_dot_bar_phi)' + (R_s*A_bar'*inv(R_x)*A_bar*R_s).*conj(A_dot_bar_theta'*inv(R_x)*A_dot_bar_phi)');
    
    %% 构造F_theta_phi
    F_theta_phi = 2*L*real((R_s*(A_bar')*inv(R_x)*A_dot_bar_phi).*conj(R_s*(A_bar')*inv(R_x)*A_dot_bar_theta)' + (R_s*A_bar'*inv(R_x)*A_bar*R_s).*conj(A_dot_bar_phi'*inv(R_x)*A_dot_bar_theta)');
    
    %% 构造F_a1_a1
    F_a1_a1 = 2*L*real(trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'));
    
    %% 构造F_a1_a2
    F_a1_a2 = 2*L*real(trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'));
    
    %% 构造F_a1_a3
    F_a1_a3 = 2*L*real(trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'));
    
    %% 构造F_a2_a1
    F_a2_a1 = 2*L*real(trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'));
    
    %% 构造F_a2_a2
    F_a2_a2 = 2*L*real(trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'));
    
    %% 构造F_a2_a3
    F_a2_a3 = 2*L*real(trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'));
    
    %% 构造F_a3_a1
    F_a3_a1 = 2*L*real(trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'));
    
    %% 构造F_a3_a2
    F_a3_a2 = 2*L*real(trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'));
    
    %% 构造F_a3_a3
    F_a3_a3 = 2*L*real(trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'));
    
    %% 构造F_b1_b1
    F_b1_b1 = 2*L*real(trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'));
    
    %% 构造F_b1_b2
    F_b1_b2 = 2*L*real(trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'));
    
    %% 构造F_b1_b3
    F_b1_b3 = 2*L*real(trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'));
    
    %% 构造F_b2_b1
    F_b2_b1 = 2*L*real(trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'));
    
    %% 构造F_b2_b2
    F_b2_b2 = 2*L*real(trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'));
    
    %% 构造F_b2_b3
    F_b2_b3 = 2*L*real(trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'));
    
    %% 构造F_b3_b1
    F_b3_b1 = 2*L*real(trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'));
    
    %% 构造F_b3_b2
    F_b3_b2 = 2*L*real(trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'));
    
    %% 构造F_b3_b3
    F_b3_b3 = 2*L*real(trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'));
    
    %% 构造F_a1_b1
    F_a1_b1 = 2*L*real(trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'));
    
    %% 构造F_a1_b2
    F_a1_b2 = 2*L*real(trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'));
    
    %% 构造F_a1_b3
    F_a1_b3 = 2*L*real(trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'));
    
    %% 构造F_a2_b1
    F_a2_b1 = 2*L*real(trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'));
    
    %% 构造F_a2_b2
    F_a2_b2 = 2*L*real(trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'));
    
    %% 构造F_a2_b3
    F_a2_b3 = 2*L*real(trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'));
    
    %% 构造F_a3_b1
    F_a3_b1 = 2*L*real(trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'));
    
    %% 构造F_a3_b2
    F_a3_b2 = 2*L*real(trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'));
    
    %% 构造F_a3_b3
    F_a3_b3 = 2*L*real(trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'));
    
    %% 构造F_b1_a1
    F_b1_a1 = 2*L*real(trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'));
    
    %% 构造F_b1_a2
    F_b1_a2 = 2*L*real(trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'));
    
    %% 构造F_b1_a3
    F_b1_a3 = 2*L*real(trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'));
    
    %% 构造F_b2_a1
    F_b2_a1 = 2*L*real(trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'));
    
    %% 构造F_b2_a2
    F_b2_a2 = 2*L*real(trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'));
    
    %% 构造F_b2_a3
    F_b2_a3 = 2*L*real(trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'));
    
    %% 构造F_b3_a1
    F_b3_a1 = 2*L*real(trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'));
    
    %% 构造F_b3_a2
    F_b3_a2 = 2*L*real(trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'));
    
    %% 构造F_b3_a3
    F_b3_a3 = 2*L*real(trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar')+trace(inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'));
    
    %% 构造F_theta_a1
    F_theta_a1 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'*inv(R_x)*A_dot_bar_theta));
    
    %% 构造F_theta_a2
    F_theta_a2 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'*inv(R_x)*A_dot_bar_theta));
    
    %% 构造F_theta_a3
    F_theta_a3 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'*inv(R_x)*A_dot_bar_theta));
    
    %% 构造F_theta_b1
    F_theta_b1 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'*inv(R_x)*A_dot_bar_theta));
    
    %% 构造F_theta_b2
    F_theta_b2 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'*inv(R_x)*A_dot_bar_theta));
    
    %% 构造F_theta_b3
    F_theta_b3 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'*inv(R_x)*A_dot_bar_theta));
    
    %% 构造F_phi_a1
    F_phi_a1 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a1'*inv(R_x)*A_dot_bar_phi));
    
    %% 构造F_phi_a2
    F_phi_a2 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a2'*inv(R_x)*A_dot_bar_phi));
    
    %% 构造F_phi_a3
    F_phi_a3 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_a3'*inv(R_x)*A_dot_bar_phi));
    
    %% 构造F_phi_b1
    F_phi_b1 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b1'*inv(R_x)*A_dot_bar_phi));
    
    %% 构造F_phi_b2
    F_phi_b2 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b2'*inv(R_x)*A_dot_bar_phi));
    
    %% 构造F_phi_b3
    F_phi_b3 = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(R_s*A_bar'*inv(R_x)*A_bar*R_s*A_dot_bar_b3'*inv(R_x)*A_dot_bar_phi));
    
    %% 构造F_a1_theta
    F_a1_theta = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(A_dot_bar_theta'*inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_a2_theta
    F_a2_theta = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(A_dot_bar_theta'*inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_a3_theta
    F_a3_theta = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(A_dot_bar_theta'*inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_b1_theta
    F_b1_theta = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(A_dot_bar_theta'*inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_b2_theta
    F_b2_theta = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(A_dot_bar_theta'*inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_b3_theta
    F_b3_theta = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_theta)+diag(A_dot_bar_theta'*inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_a1_phi
    F_a1_phi = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(A_dot_bar_phi'*inv(R_x)*A_dot_bar_a1*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_a2_phi
    F_a2_phi = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(A_dot_bar_phi'*inv(R_x)*A_dot_bar_a2*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_a3_phi
    F_a3_phi = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(A_dot_bar_phi'*inv(R_x)*A_dot_bar_a3*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_b1_phi
    F_b1_phi = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(A_dot_bar_phi'*inv(R_x)*A_dot_bar_b1*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_b2_phi
    F_b2_phi = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(A_dot_bar_phi'*inv(R_x)*A_dot_bar_b2*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 构造F_b3_phi
    F_b3_phi = 2*L*real(diag(R_s*A_bar'*inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_dot_bar_phi)+diag(A_dot_bar_phi'*inv(R_x)*A_dot_bar_b3*R_s*A_bar'*inv(R_x)*A_bar*R_s))';
    
    %% 汇总求解FIM
    F_All = [F_theta_theta,F_theta_phi,F_theta_a1,F_theta_a2,F_theta_a3,F_theta_b1,F_theta_b2,F_theta_b3;
             F_phi_theta,F_phi_phi,F_phi_a1,F_phi_a2,F_phi_a3,F_phi_b1,F_phi_b2,F_phi_b3;
             F_a1_theta,F_a1_phi,F_a1_a1,F_a1_a2,F_a1_a3,F_a1_b1,F_a1_b2,F_a1_b3;
             F_a2_theta,F_a2_phi,F_a2_a1,F_a2_a2,F_a2_a3,F_a2_b1,F_a2_b2,F_a2_b3;
             F_a3_theta,F_a3_phi,F_a3_a1,F_a3_a2,F_a3_a3,F_a3_b1,F_a3_b2,F_a3_b3;
             F_b1_theta,F_b1_phi,F_b1_a1,F_b1_a2,F_b1_a3,F_b1_b1,F_b1_b2,F_b1_b3;
             F_b2_theta,F_b2_phi,F_b2_a1,F_b2_a2,F_b2_a3,F_b2_b1,F_b2_b2,F_b2_b3;
             F_b3_theta,F_b3_phi,F_b3_a1,F_b3_a2,F_b3_a3,F_b3_b1,F_b3_b2,F_b3_b3];
    
    G = inv(F_All);
    G_diag = diag(G);
    
    %% 计算DOA的CRB
    temp_CRB = 0;
    theta_CRB = 0;
    phi_CRB = 0;
    for idtemp = 1:2*K
        temp_CRB = temp_CRB + G_diag(idtemp);
    end
    for idtemp = 1:K
        theta_CRB = theta_CRB + G_diag(idtemp);
    end
    for idtemp = K+1:2*K
        phi_CRB = phi_CRB + G_diag(idtemp);
    end
    CRB_DOA(snr_index) = sqrt(temp_CRB/(2*K)) * 180/pi; % 转换为度
    CRB_theta(snr_index) = sqrt(theta_CRB/(K)) * 180/pi; % 转换为度
    CRB_phi(snr_index) = sqrt(phi_CRB/(K)) * 180/pi; % 转换为度
    %% 计算互耦参数的CRB
    temp_CRB_C = 0;
    for idtemp = 2*K+1:length(G_diag)
        temp_CRB_C = temp_CRB_C + G_diag(idtemp);
    end
    CRB_C(snr_index) = sqrt(temp_CRB_C)/norm_C;
end

%% 绘图
figure

% DOA的CRB vs SNR
semilogy(SNR_list, CRB_DOA, '-o', 'LineWidth', 1.5);hold on
semilogy(SNR_list, CRB_phi, '-*', 'LineWidth', 1.5);hold on
semilogy(SNR_list, CRB_theta, '-^', 'LineWidth', 1.5)
grid on
xlabel('SNR (dB)', 'FontSize', 12)
ylabel('CRB for DOA (degrees)', 'FontSize', 12)
title('DOA估计的CRB界 vs SNR', 'FontSize', 14)
set(gca, 'FontSize', 11)
ylim([1e-2,1])
% 互耦参数的CRB vs SNR
figure
semilogy(SNR_list, CRB_C, '-s', 'LineWidth', 1.5)
grid on
xlabel('SNR (dB)', 'FontSize', 12)
ylabel('CRB for Coupling Parameters (normalized)', 'FontSize', 12)
title('互耦参数估计的CRB界 vs SNR', 'FontSize', 14)
set(gca, 'FontSize', 11)

%% 保存结果
filename = sprintf('CRB_vs_SNR_L=%d_N=%dx%d_K=%d.mat', L, N_x, N_y, K);
save(filename, 'CRB_DOA', 'CRB_C', 'SNR_list', 'L', 'N_x', 'N_y', 'K');

toc