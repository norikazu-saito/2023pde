% Error observation
% the implicit theta FDM forn the heat equation
% u_t=ku_xx+f with 0-DBC
% h = 1 / (Nx+1)  
% tau = Tmax/Nt
function[tau,error1,error2] = heat2(Nx,Nt,theta,Tmax,givenfunc,id)
%%% if id >= 0 then plot graph
%%% if id < 0  then no graph
% diffusion coef
k=1;
%% space interval and space mesh
N=Nx;
a = 0; b = 1; h = (b - a) / (N+1);
x = (a+h : h: b-h).'; xx = (a: h: b).'; 
%% boundary values
ua=0; ub=0;
%% time increment  
tau = Tmax/Nt; 
lam = k*tau/(h^2); nmax=floor(Tmax/tau); 
number = 40; step = floor(max(1, nmax/number));
%% discrete Laplacian
A = 2 * eye(N, N) - diag(ones(N - 1, 1), 1) - diag(ones(N - 1, 1), -1); 
%% theta-method
H = sparse(eye(N, N) + theta * lam * A);
K = sparse(eye(N, N) - (1 - theta) * lam * A);
%% LU factorization
[L, U, P] = lu(H);
%% initial time 
tnow = 0.0; 
%% set initial value
u = givenfunc(0,x,0); uu=[ua;u;ub]; 
%% draw
if id>=0 
    %% for 2d draw
    figure(1); hold on; plot(xx,uu,'c');
    %% for 3d draw
    figure(2); hold on; tsp=tnow*ones(1,N+2); plot3(xx,tsp,uu,'c');
end
%% iteration 
time=[]; error1 = -1; error2 = -1;  
for n=1:nmax
    tpast=tnow; tnow = n*tau; 
    f = K*u + tau*((1-theta)*givenfunc(1, x, tpast) + theta*givenfunc(1, x, tnow));
    u = (U \ (L \ (P * f)));
    uu=[ua;u;ub]; 
    % error
    err = u - givenfunc(2, x, tnow); 
    error1 = max(error1, norm(err,inf));
    error2 = max(error2, norm(err,2)); error2=error2*sqrt(h);
    % draw 
    if id>=0 
        if rem(n, step)==0
        % for 2d & 3d draw 
        figure(1); plot(xx,uu,'c');
        figure(2); tsp=tnow*ones(1,N+2); plot3(xx,tsp,uu,'c');
        end
    end
end
if id>=0 
% decoration of figure windows
    figure(1);xlabel('x');ylabel('u');grid on;saveas(1,'heat2a.pdf');
    figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
    view(60,15);saveas(2,'heat2b.pdf') 
end
%%% end of program
end