function [u_hat]= fun_1D_search_u(a_get,u_grid,Posd,lambda)
cu=u_grid(2)-u_grid(1);
f_fun_cu=zeros(length(u_grid),1);
for i=1:length(u_grid)
    a_grid_cu=exp(-1j*2*pi/lambda*Posd*(u_grid(i)));
    f_fun_cu(i)=abs(a_get'*a_grid_cu)/(norm(a_get).^2)/(norm(a_grid_cu).^2);
end
[~,idxcu]=max(f_fun_cu);
theta_cu=u_grid(idxcu(1));
xi=cu/100;
grid_xi=theta_cu-cu:xi:theta_cu+cu;
f_fun_xi=zeros(length(grid_xi),1);
for i=1:length(grid_xi)
    a_grid_xi=exp(-1j*2*pi/lambda*Posd*(grid_xi(i)));
    f_fun_xi(i)=abs(a_get'*a_grid_xi)/(norm(a_get).^2)/(norm(a_grid_xi).^2);
end
[~,idxcu]=max(f_fun_xi);
u_hat=grid_xi(idxcu(1));