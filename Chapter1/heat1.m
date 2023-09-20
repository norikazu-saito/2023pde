% The explicit FDM for the heat equation
% u_t=ku_xx with 0-DBC
function heat1(N, lam, Tmax, givenfunc)
% diffusion coef
k=1;
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
u = givenfunc(0,x,0); uu=[ua;u;ub]; 
%% draw
%% for 2d draw
figure(1); hold on; plot(xx,uu,'r');
%% for 3d draw 
figure(2); hold on; tsp=tnow*ones(1,N+2); plot3(xx,tsp,uu,'r');
%% iteration 
time=[];
for n=1:nmax   
    tpast = tnow; tnow = n*tau; u = K*u + tau*givenfunc(1, x, tpast); uu=[ua;u;ub]; 
    %% for 2d & 3d draw
    if rem(n, step)==0
        figure(1); plot(xx,uu,'b');
        figure(2); tsp=tnow*ones(1,N+2); plot3(xx,tsp,uu,'b');
    end
end
% decoration of figure windows
figure(1);xlabel('x');ylabel('u');grid on;saveas(1,'heat1a.pdf');
figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
view(60,15);saveas(2,'heat1b.pdf'); 
%%% end of program
end
