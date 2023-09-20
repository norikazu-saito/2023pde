%% Explicit scheme for Cahn-Hilliard equation
function ch01(p, q, r, N, tau, Tmax)
%% space interval and space mesh
a = 0; b = 1; h = (b - a) / N; x = (a: h: b).'; 
%% CFL number and time increment
nmax=floor(Tmax/tau); number = 30; step = floor(max(1, nmax/number));
%% discrete operator
A = 2 * eye(N+1, N+1) - diag(ones(N, 1), 1) - diag(ones(N, 1), -1); 
A(1,2)=-2; A(N+1,N)=-2; A = (1/(h^2))*sparse(A); 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set initial value
tnow = 0.0; u = ch_init1(x);
%% draw
figure(1); hold on; plot(x,u,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iteration 
n=0; kmax=30; 
while (tnow < Tmax) 
    n = n + 1; tnow = n*tau; 
    w = u - tau*A*(p*u + r*u.^3 - q*A*u);
    u = w;  
    %%%%%%%%%%%%%%%%%%%%%
    % draw 
    figure(1); plot(x,u,'b','LineWidth',2);
end
% decoration of figure windows
figure(1);xlabel('x');ylabel('u');grid on;%pbaspect([2 1 1]);
saveas(1,'ch01a.pdf');
end
%%%%
%%% initial value
function w = ch_init1(x)
   w = 0.1*sin(2*pi*x) + 0.01*cos(4*pi*x) + 0.06*sin(4*pi*x) + 0.02*cos(10*pi*x); 
end
%%%