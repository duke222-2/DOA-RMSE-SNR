clc;clear;close all;
% randn('seed', 0);
% rand('seed', 0);
K=2;
M=10;
T=200;
c=3e8;
f=3e9;
lambda=c/f;
d=lambda/2;
Pos=d*(0:M-1);
interval =0:1:10;
angle=zeros(length(interval),1);
grid_theta=0:0.1:180;
SNR=20;
numTrials=300;
theta_MUSIC=zeros(numTrials,K);
theta_ESPRIT=zeros(numTrials,K);
%%%%%%%%%%%%%%%%%%%%
A=zeros(M,K);
recPb_array_MUSIC=zeros(length(interval),1);
recPb_array_ESPRIT=zeros(length(interval),1);
%%%%%%%%%%%%%%%%%%%%
h = waitbar(0, '开始蒙特卡洛仿真...');
for idx=1:length(interval)
    angle=interval(idx);
    theta_true=[30,30+angle];
    for k=1:K
        A(:,k)=exp(1j*2*pi/lambda*Pos.'*cosd(theta_true(k)));
    end
    for trial=1:numTrials
        waitbar((idx-1)/length(interval) + (trial/numTrials)*(1/length(interval)), h, ...
            sprintf('SNR = %.1f dB | 试验 %d/%d | 总体进度: %.1f%%', ...
            SNR, trial, numTrials, 100*((idx-1)/length(interval) + (trial/numTrials)*(1/length(interval)))));
        S=randn(K,T)+1j*randn(K,T);
        Y=A*S;
        Y=awgn(Y,SNR,"measured");
        theta_ESPRIT(trial,:)=fun_esprit_1D(Y,K,M,T,lambda,d);
        xi_len=10;
        theta_MUSIC(trial,:)=fun_MUSIC_1D(Y,K,Pos,lambda,grid_theta,xi_len);
    end
    recPb_array_ESPRIT(idx) = func_recPb(theta_true,theta_ESPRIT,angle);
    recPb_array_MUSIC(idx) = func_recPb(theta_true,theta_MUSIC,angle);
end
%%%%%%%%%%%%%%%%%%%%
base_name = 'resolution_1D_DOA';
filename = strcat(base_name, ...
    '_K',     num2str(K), ...
    '_M',     num2str(M), ...
    '_T',     num2str(T), ...
    '_Trials',num2str(numTrials), ...
    '_angle',   num2str(min(interval)), 'to', num2str(max(interval)), '.mat');
save(filename, 'recPb_array_MUSIC','recPb_array_ESPRIT','interval');
%%%%%%%%%%%%%%%%%%%%
set(0,'defaultAxesFontSize',13);
set(0,'defaultTextFontSize',13);
set(0,'defaultLegendFontSize',13);
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');
%%%%%%%%%%%%%%%%%%%%
default_colors=lines(2);

close(h)

figure

plot(interval,recPb_array_ESPRIT, '-o', 'LineWidth', 1.5,...
    'Color', default_colors(1, :));hold on;

plot(interval,recPb_array_MUSIC, '-*', 'LineWidth', 1.5,...
    'Color', default_colors(2, :));hold on;

xlim([min(interval),max(interval)]);

grid on
title(sprintf('MUSIC&ESPRIT'),"FontSize",13);
ylabel('Resolution');
xlabel('\delta')
legend('ESPRIT','MUSIC','Location', 'northeast');

exportgraphics(gcf,'resolution_1D_DOA.pdf','ContentType','vector');