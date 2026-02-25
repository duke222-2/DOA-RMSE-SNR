function doa=fun_esprit_1D(Y,K,M,T,lambda,d)
%输出是行向量
R=Y*Y'/T;%M*M
[E,D]=eig(R);
[~,I]=sort(diag(D),'descend');
E=E(:,I);
E=E(:,1:K);%M*K
E1=E(1:M-1,:);
E2=E(2:M,:);

Z=[E1,E2];%(M-1)*2K
R_tls=Z'*Z;%2K*2K
[E_t,D_t]=eig(R_tls);
[~,I_t]=sort(diag(D_t));
E_t=E_t(:,I_t);
F=E_t(:,1:K);%2K*K
F0=F(1:K,:);
F1=F(K+1:2*K,:);
Psi=-F0/F1;
z=eig(Psi);
doa=acos(angle(z)/(2*pi)*lambda/d)/pi*180;% d*sin(theta)=φ/(2π)*lambda
doa=sort(real(doa));
doa=doa.';
end

