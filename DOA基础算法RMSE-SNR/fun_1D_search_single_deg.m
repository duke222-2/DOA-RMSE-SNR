function [theta_hat]= fun_1D_search_single_deg(a_get,theta_grid,Posd,lambda)
cu=theta_grid(2)-theta_grid(1);
f_fun_cu=zeros(length(theta_grid),1);
for i=1:length(theta_grid)
    a_grid_cu=exp(1j*2*pi/lambda*Posd*cosd(theta_grid(i)));
    f_fun_cu(i)=abs(a_get'*a_grid_cu)/(norm(a_get).^2)/(norm(a_grid_cu).^2);
end
[~,idxcu]=max(f_fun_cu);
theta_cu=theta_grid(idxcu(1));
xi=cu/100;
grid_xi=theta_cu-cu:xi:theta_cu+cu;
f_fun_xi=zeros(length(grid_xi),1);
for i=1:length(grid_xi)
    a_grid_xi=exp(1j*2*pi/lambda*Posd*cosd(grid_xi(i)));
    f_fun_xi(i)=abs(a_get'*a_grid_xi)/(norm(a_get).^2)/(norm(a_grid_xi).^2);
end
[~,idxcu]=max(f_fun_xi);
theta_hat=grid_xi(idxcu(1));