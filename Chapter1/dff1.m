function dff1(N, lam, r, Tmax)
% if r=0, the scheme solves the wave equation. 
% if r=1, the scheme solves the heat equation.
%% parameters 
a=0; b=1; h=(b - a)/(N+1); x=(a+h:h:b-h).'; xx=(a:h:b).'; ua=0; ub=0; 
%% time increment  
tau=lam*h^(r+1); nmax=floor(Tmax/tau); number=60; step=floor(max(1,nmax/number)); 
p= 2*lam/(1+2*lam); q= (1-2*lam)/(1+2*lam);  
%% 
B = sparse(diag(ones(N-1,1),1)+diag(ones(N-1,1),-1));
%% set initial value 
figure(2); hold on; 
tnow=0; v=exact(x,tnow); 
uplt=[ua;v;ub]; tsp=tnow*ones(1,N+2); plot3(xx,tsp,uplt,'r');
tnow=tau; u=exact(x,tau); 
uplt=[ua;u;ub]; tsp=tnow*ones(1,N+2); plot3(xx,tsp,uplt,'b');
%% iteration
for n=2:nmax   
    unew = p*B*u+q*v; uplt=[ua;unew;ub]; v=u; u=unew; tnow = n*tau; 
    %% for 3d draw
    if rem(n, step)==0
        figure(2); tsp=tnow*ones(1,N+2); 
        %plot3(xx,tsp,uplt,'c');
        plot3(xx,tsp,uplt,'b');
    end
end
% decoration
figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
view(45,45);saveas(2,'dff1.pdf'); 
end
%%%
function z=exact(x,t)
    z=exp(-pi^2*t).*sin(pi*x); 
    %z=sin(4*pi*x)*(1+t);
end
