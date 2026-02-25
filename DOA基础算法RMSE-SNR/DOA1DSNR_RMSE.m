clc;clear;close all;
randn('seed', 0);
rand('seed', 0);
K=3;
M=16;
T=100;
c=3e8;
f=3e9;
lambda=c/f;
d=lambda/2;
Pos=d*(0:M-1);
theta=[30,60,135];
grid_theta=0:1:180;
theta_MUSIC=zeros(1,K);
theta_ESPRIT=zeros(1,K);
SNR_list=-0:5:20;
numTrials=200;
mse_MUSIC = zeros(numTrials, length(SNR_list));
mse_ESPRIT = zeros(numTrials, length(SNR_list));
A=zeros(M,K);
for k=1:K
    A(:,k)=exp(1j*2*pi/lambda*Pos.'*cosd(theta(k)));
end
h = waitbar(0, '开始蒙特卡洛仿真...');
for idx=1:length(SNR_list)
    SNR=SNR_list(idx);
    for trial=1:numTrials
        waitbar((idx-1)/length(SNR_list) + (trial/numTrials)*(1/length(SNR_list)), h, ...
            sprintf('SNR = %.1f dB | 试验 %d/%d | 总体进度: %.1f%%', ...
            SNR, trial, numTrials, 100*((idx-1)/length(SNR_list) + (trial/numTrials)*(1/length(SNR_list)))));
                S=randn(K,T)+1j*randn(K,T);
        Y=A*S;
        Y=awgn(Y,SNR,"measured");
        xi_len=100;
        theta_MUSIC=fun_MUSIC_1D(Y,K,Pos,lambda,grid_theta,xi_len);
        theta_ESPRIT=fun_esprit_1D(Y,K,M,T,lambda,d);
        mse_MUSIC(trial,idx)=mean(abs(sort(theta) - sort(theta_MUSIC)).^2);
        mse_ESPRIT(trial,idx)=mean(abs(sort(theta) - sort(theta_ESPRIT)).^2);
    end
end
RMSE_MUSIC = sqrt(mean(mse_MUSIC, 1));
RMSE_ESPRIT=sqrt(mean(mse_ESPRIT,1));
base_name = 'MUSICvsESPRIT';

filename = strcat(base_name, ...
    '_K',     num2str(K), ...
    '_M',     num2str(M), ...
    '_T',     num2str(T), ...
    '_Trials',num2str(numTrials), ...
    '_SNR',   num2str(min(SNR_list)), 'to', num2str(max(SNR_list)), '.mat');
save(filename, 'RMSE_MUSIC','RMSE_ESPRIT','SNR_list');


%%
set(0,'defaultAxesFontSize',13);
set(0,'defaultTextFontSize',13);
set(0,'defaultLegendFontSize',13);
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');

default_colors=lines(2);

close(h)

figure

semilogy(SNR_list,RMSE_ESPRIT, '-o', 'LineWidth', 1.5,...
    'Color', default_colors(1, :));hold on;
semilogy(SNR_list,RMSE_MUSIC, '-*', 'LineWidth', 1.5,...
    'Color', default_colors(2, :));hold on;

ylim([1e-2,1e2]);
xlim([min(SNR_list),max(SNR_list)]);

grid on
title(sprintf('MUSIC&ESPRIT'),"FontSize",13);
ylabel('RMSE (°)');
xlabel('SNR (dB)');
legend('ESPRIT','MUSIC',  'Location', 'northeast');

exportgraphics(gcf,'MUSIC&ESPRIT_Performance.pdf','ContentType','vector');