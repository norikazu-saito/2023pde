%% Discrete variational derivative method for Cahn-Hilliard equation
function ch11(p, q, r, N, tau, Tmax, tol, id)
%%% if id >= 0 then plot graph
%%% if id < 0  then no graph
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
if id>=0 
    %% for 2d draw
    figure(1); hold on; plot(x,u,'r');
    %% for 3d draw
    figure(2); hold on; tsp=tnow*ones(1,N+1); plot3(x,tsp,u,'r');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iteration 
n=0; tau1=1.0/tau; kmax=30; 
while (tnow < Tmax) 
    n = n + 1; tnow = n*tau; 
    %%% Newton iteration 
    k = 1; dif = tol + 1; 
    w = u; mag=norm(u, inf);
    while (k < kmax) && (dif > tol)
        ww = 0.5*(w+u); ww2 = 0.5*(w.^2+u.^2);
        F = tau1*(w-u) + A*(p*ww + r*ww.*ww2 - q*A*ww); 
        ww3 = 0.5*ww2 + 0.5*(w+u).*w; 
        DF = tau1*eye(N+1,N+1) + 0.5*p*A + r*A*diag(ww3) -0.5*q*A*A;
        [L, U, P] = lu(DF);
        v = (U \ (L \ (P * F)));
        w = w - v;
        dif = norm(v, inf)/mag;
        k = k + 1;
    end
    % Newton iteration converge??
    if k >= kmax
        disp("Newton method does not converge!"); break;
    end 
    u = w;  
    %%%%%%%%%%%%%%%%%%%%%
    % draw 
    if id>=0 
        if rem(n, step)==0
        % for 2d & 3d draw 
        figure(1); plot(x,u,'b');
        figure(2); tsp=tnow*ones(1,N+1); plot3(x,tsp,u,'b');
        end
    end
end
if id>=0 
% decoration of figure windows
    figure(1);xlabel('x');ylabel('u');grid on;saveas(1,'ch11a.pdf');
    figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
    view(25,54);saveas(2,'ch11b.pdf') 
end
end
%%%%
%%% initial value
function w = ch_init1(x)
   w = 0.1*sin(2*pi*x) + 0.01*cos(4*pi*x) + 0.06*sin(4*pi*x) + 0.02*cos(10*pi*x); 
end
%%%