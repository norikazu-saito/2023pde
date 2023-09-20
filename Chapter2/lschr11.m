function [conservation]=lschr11(N, theta, Tmax)
% parameters
a=0; b=2*pi; h=(b-a)/(N+1); x=(a+h:h:b-h).'; xx=(a:h:b).';
ua=0+1i*0; ub=0+1i*0; tau=0.1*h; lambda = tau/h/h; nmax=floor(Tmax/tau); 
number = 40; step = floor(max(1, nmax/number)); conservation=[];
% discrete Laplacian
Phi=diag(func_phi(x)); 
MM=2*eye(N, N)-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1);
H=sparse(eye(N, N)+1i*theta*lambda*MM+1i*theta*tau*Phi); 
K=sparse(eye(N, N)-1i*(1-theta)*lambda*MM-1i*(1-theta)*tau*Phi); 
MM=(1/(h^2))*sparse(MM); 
% LU factorization
[L,U,P]=lu(H);
% set initial value
tnow=0; u=func_a(x); uu=[ua;u;ub]; 
mass=h*norm(u,2)^2; energy=(h/2)*(dot(MM*u,u)+dot(Phi*u,u)); 
conservation=[conservation;tnow,mass,energy]; 
figure(1); hold on; tsp=tnow*ones(1,N+2); plot3(xx,tsp,abs(uu),'r');
% iteration 
for n=1:nmax
    tnow=n*tau; f=K*u; w=(U\(L\(P * f))); u=w;
    mass=h*norm(u,2)^2; energy=(h/2)*(dot(MM*u,u)+dot(Phi*u,u)); 
    conservation=[conservation;tnow,mass,energy]; 
    if rem(n, step)==0
        figure(1); tsp=tnow*ones(1,N+2); uu=[ua;u;ub]; plot3(xx,tsp,abs(uu),'b');
    end
end
conservation=real(conservation); 
% decoration
figure(1);xlabel('x');ylabel('t');zlabel('abs(u)');grid on; 
view(60,52);saveas(1,'lschr11.pdf') 
end
% initial value
function w = func_a(x)
    w = sin(x).^2 + 1i * sin(x).^2;
    %w = cos(x).^2 + 1i * cos(x).^2;
end
% potential 
function w = func_phi(x)
    w = x.^2 + 1; 
end
