function theta_est=fun_MUSIC_1D(Y,K,Pos,lambda,grid_theta,xi_len)
T=size(Y,2);
R      =  Y*Y'/T;
[EV,D] =  eig(R);
[~,I] = sort(diag(D),'descend');
EV= EV(:,I);
EN=EV(:,K+1:end);

cu=grid_theta(2)-grid_theta(1);
P_MUSIC = zeros(1, length(grid_theta));
for ii = 1:length(grid_theta)
    a = exp(1j*2*pi/lambda*cosd(grid_theta(ii))*Pos.');
    P_MUSIC(ii) =  1/abs(a'*(EN*EN')*a);
end
[peaks, loc] = findpeaks(P_MUSIC, grid_theta);

if length(peaks)<K
    loc(end+1:K)=0;
else
    [~,idx]=sort(peaks,"descend");
    loc=loc(idx(1:K));
end

xi=cu/xi_len;
fine_angles = zeros(1,K);
for ii=1:K
    if loc(ii)==0
        fine_angles(ii)=0;
        continue;
    end
    regions=loc(ii)-cu:xi:loc(ii)+cu;
    P_fine=zeros(1,length(regions));
    for jj=1:length(regions)
        a = exp(1j*2*pi/lambda*cosd(regions(jj))*Pos.');
        P_fine(jj) = 1 / (a' * (EN * EN') * a);
    end
    [~,idx_fine]=max(P_fine);
    fine_angles(ii)=regions(idx_fine);
end
theta_est = sort(fine_angles);
end
