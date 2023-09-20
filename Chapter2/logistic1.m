% The explicit FDM for a semilinear equation
% u_t=ku_xx + ep(1-u)u with 0-DBC
function [r,rho] = logistic1(k, ep, N, lam, Tmax,func_iv)
%% space interval and space mesh
a = 0; b = 1; h = (b - a) / (N+1);
x = (a+h : h: b-h).'; xx = (a: h: b).'; 
%% boundary values
ua=0; ub=0;
%% time increment  
tau = lam*h^2/k; nmax=floor(Tmax/tau); 
number = 40; step = floor(max(1, nmax/number)); 
%% discrete Laplacian
A = 2 * eye(N, N) - diag(ones(N - 1, 1), 1) - diag(ones(N - 1, 1), -1); 
K = sparse(eye(N, N) - lam * A);
%% initial time 
tnow = 0.0; 
%% set initial value
u = func_iv(x); uu=[ua;u;ub]; 
%% draw
%% for 2d draw
figure(1); hold on; plot(xx,uu,'r');
%% for 3d draw 
figure(2); hold on; tsp=tnow*ones(1,N+2); plot3(xx,tsp,uu,'r');
%% iteration 
time=[];
for n=1:nmax   
    tpast = tnow; tnow = n*tau; u = K*u + tau*ep*(1-u).*u; uu=[ua;u;ub]; 
    %% for 2d & 3d draw
    if rem(n, step)==0
        figure(1); plot(xx,uu,'b');
        figure(2); tsp=tnow*ones(1,N+2); plot3(xx,tsp,uu,'b');
    end
end
% decoration of figure windows
figure(1);xlabel('x');ylabel('u');grid on;saveas(1,'logistic1a.pdf');
figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
view(33,44);saveas(2,'logistic1b.pdf'); 
%%% end of program
r = ep / (k*pi^2); rho = tau*(ep + 2*k/(h^2));
end
