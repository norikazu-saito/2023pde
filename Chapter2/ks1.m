function ks1(N,Tmax,mass,id)
% parameters 
a=0; b=1;xx=linspace(a,b,N+1)'; x=arrange1(xx);  
h=xx(2)-xx(1); lam=0.5; gamma=1; alpha=1; mu=10; ep=0.9;  
%%
A=2*eye(N+1,N+1)-diag(ones(N,1),1)-diag(ones(N,1),-1);
A(1,2)=-2;A(N+1,N)=-2; A=sparse((1/h^2)*A+gamma*eye(N+1,N+1)); 
[L,U,P]=lu(A);
%% time increment  
tau0=lam*h^2; nmax=floor(Tmax/tau0); number=40; step=floor(max(1,nmax/number)); 
%% set initial value
tnow = 0.0; u = init(id,mass,x); 
figure(2); hold on; tsp=tnow*ones(1,N); plot3(x,tsp,u,'r');
%% iteration 
for n=1:nmax
    f=alpha*arrange2(u); v=(U\(L\(P*f))); b=(mu/h)*arrange3(v); 
    q=arrange4(b,u,h); tau=min(tau0, ep*h^2/(2+h*norm(b,inf))); 
    unew=u-(tau/h)*arrange3(q); u=unew; tnow = n*tau; 
    if rem(n, step)==0
        figure(2); tsp=tnow*ones(1,N); plot3(x,tsp,u,'b');
    end
end
figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
view(30,60);saveas(2,'ks1.pdf');
end
%%%%
% initial value
function w=init(id, m, x)
del=0.1; 
switch id
    case 1  
        I=find(x>0.25 & x <=0.45); J=find(x>0.55 & x <=0.75); 
        z=zeros(size(x)); z(I)=2.5; z(J)=2.5;
        w=(m-del)*z+del; 
    case 2 
        I=find(x>0.1 & x <=0.3); J=find(x>0.7 & x <=0.9); 
        z=zeros(size(x)); z(I)=2.5; z(J)=2.5;
        w=(m-del)*z+del; 
    case 3 
        z = zeros(size(x)); 
end
end
function x=arrange1(xx)
    x=0.5*(xx+circshift(xx,-1)); x=x(1:end-1);
end 
function f=arrange2(u)  
    f=0.5*(u+circshift(u,-1)); f(end)=u(end); f=[u(1);f];
end 
function b=arrange3(v)
    b=v-circshift(v,1); b=b(2:size(v));
end 
function qq=arrange4(b,u,h)
    bp=max(0,b);bm=max(0,-b); uu=circshift(u,-1); bbm=circshift(bm,-1); 
    qq=-(1/h*ones(size(b))+bbm).*uu+(1/h*ones(size(b))+bp).*u; 
    qq=[0;qq]; qq(end)=0; 
end 
