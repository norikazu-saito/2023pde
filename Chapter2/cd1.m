function cd1(N,Tmax)
% parameters 
a=0; b=1; xx=linspace(a,b,N+1)'; x=xx(1:end-1); h=xx(2)-xx(1); 
lam=0.4; tau=lam*h*h; 
%%
A=2*eye(N,N)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);A(1,N)=-1;A(N,1)=-1;
B=diag(ones(N-1,1),1)-diag(ones(N-1,1),-1); B(1,N)=-1;B(N,1)=1;
A=sparse(eye(N,N)-lam*A); B=(tau/(2*h))*sparse(B);  
%% time increment  
tau=lam*h^2; nmax=floor(Tmax/tau); number=40; step=floor(max(1,nmax/number)); 
%% set initial value
tnow = 0.0; u = init(x); uu=[u;u(end)]; 
figure(2); hold on; tsp=tnow*ones(1,N+1); plot3(xx,tsp,uu,'r');
%% iteration 
for n=1:nmax 
    unew=A*u-B*(conv(x,tnow).*u); u=unew; tnow = n*tau;
    if rem(n, step)==0
        uu=[u;u(end)];figure(2);tsp=tnow*ones(1,N+1);plot3(xx,tsp,uu,'b');
    end
end
figure(2);xlabel('x');ylabel('t');zlabel('u');grid on; 
view(30,45);saveas(2,'cd1.pdf');
end
%%%%
function w=init(x)
    w=1+sin(2*pi*x);
end
function w=conv(x,t)
    w=4*(1+cos(2*pi*x)).*(1+t).^2; 
end 