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
theta=[30,60,135];
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
        X=S.';
        RS=(X*X')/T;
        AA=kr(A_y, A_x);
        R=(AA*X)*((AA*X)')/T;
        R=R+noise_power*eye(size(R));
        %AA_theta
        A_x_thetad=zeros(size(A_x));
        A_y_thetad=zeros(size(A_y));
        AA_theta=zeros(M_x*M_y,K);
        for k=1:K
            A_x_thetad(:,k)=(1j*2*pi/lambda*Pos_x*cosd(theta(k))*sind(phi(k))).*A_x(:,k);
            A_y_thetad(:,k)=(1j*2*pi/lambda*Pos_y*cosd(theta(k))*cosd(phi(k))).*A_y(:,k);
            AA_theta(:,k)=kron(A_y_thetad(:,k),A_x(:,k)) +kron(A_y(:,k),A_x_thetad(:,k));
        end
        %AA_phi
        A_x_phid=zeros(size(A_x));
        A_y_phid=zeros(size(A_y));
        AA_phi=zeros(M_x*M_y,K);
        for k=1:K
            A_x_phid(:,k)=(1j*2*pi/lambda*Pos_x*sind(theta(k))*cosd(phi(k))).*A_x(:,k);
            A_y_phid(:,k)=(-1j*2*pi/lambda*Pos_y*sind(theta(k))*sind(phi(k))).*A_y(:,k);
            AA_phi(:,k)=kron(A_y_phid(:,k),A_x(:,k)) +kron(A_y(:,k),A_x_phid(:,k));
        end

        F_theta_theta=2*T*real((RS*AA'*inv(R)*AA_theta).*(RS*AA'*inv(R)*AA_theta).' ...
            +(RS*AA'*inv(R)*AA*RS).*(AA_theta'*inv(R)*AA_theta).');
        F_phi_phi=2*T*real((RS*AA'*inv(R)*AA_phi).*(RS*AA'*inv(R)*AA_phi).' ...
            +(RS*AA'*inv(R)*AA*RS).*(AA_phi'*inv(R)*AA_phi).');
        F_phi_theta=2*T*real((RS*AA'*inv(R)*AA_theta).*(RS*AA'*inv(R)*AA_phi).' ...
            +(RS*AA'*inv(R)*AA*RS).*(AA_theta'*inv(R)*AA_phi).');
        %FIM
        FIM=[F_theta_theta,F_phi_theta';
            F_phi_theta,F_phi_phi];
        CRB_matrix=inv(FIM);
        crb_theta(trial,idx)=sum(diag(CRB_matrix(0*K+1:1*K, 0*K+1:1*K)))/K;
        crb_phi(trial,idx)=sum(diag(CRB_matrix(1*K+1:2*K, 1*K+1:2*K)))/K;
        crb(trial,idx)=sum(diag(CRB_matrix(0*K+1:2*K,0*K+1:2*K)))/(2*K);

    end
end
CRB_theta=sqrt((mean(crb_theta,1)));
CRB_phi=sqrt((mean(crb_phi,1)));
CRB=sqrt((mean(crb,1)));

close(h)
save('CRB2D.mat','CRB_theta','CRB_phi','CRB')
figure
semilogy(SNR_list,CRB_theta*180/pi,'-o','LineWidth',1.5); hold on;
semilogy(SNR_list,CRB_phi*180/pi,'-*','LineWidth',1.5); hold on;
semilogy(SNR_list,CRB*180/pi,'-^','LineWidth',1.5); hold on;
grid on
ylim([10^(-2.7), 10^(0)]);

