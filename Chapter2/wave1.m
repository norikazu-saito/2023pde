function wave1(N, lam, Tmax)
%% parameters 
c=1; a=0; b=1; h=(b - a)/(N+1); x=(a+h:h:b-h).'; xx=(a:h:b).'; 
ua=0; ub=0; 
%% time increment  
tau=lam*h/c; nmax=floor(Tmax/tau); number=60; step=floor(max(1,nmax/number)); 
%% discrete Laplacian
A=2*eye(N,N)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1); 
B=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
K = sparse(2*eye(N, N) - lam^2 * A);
%% initial time 
tnow = 0.0; 
%% set initial value 
v = init_f(x); u=tau*init_g(x)+(1-lam^2)*v+(lam^2)/2*B*v; 
%% for 3d draw 
figure(2); hold on; 
uplt=[ua;v;ub]; tsp=tnow*ones(1,N+2); plot3(xx,tsp,uplt,'r');
%% iteration 
time=[];
for n=1:nmax   
    unew = K*u-v; uplt=[ua;unew;ub]; v=u; u=unew; tnow = n*tau; 
    %% for 3d draw
    if rem(n, step)==0
        figure(2); tsp=tnow*ones(1,N+2); plot3(xx,tsp,uplt,'b');
    end
end
% decoration
figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
view(45,45);saveas(2,'wave1.pdf'); 
end
%%%
function w=init_f(x)
    I=find(x>3/8 & x <=1/2); J=find(x>1/2 & x <=5/8); 
    w=zeros(size(x)); w(I)=8*x(I)-3; w(J)=-8*x(J)+5;
end
function w=init_g(x)
    w=zeros(size(x));
end
