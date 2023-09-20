function cd1d_uwd(m,N)
% parameters 
a = 0; b = 1; xx = linspace(a,b,N+2)'; x = xx(2:end-1); h = xx(2) - xx(1); 
%%
A = sparse(N,N); B = sparse(N,N); K = sparse(N,N); 
A = 2*eye(N,N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
B = eye(N,N) - diag(ones(N-1,1),-1);
%%
b = 2; nu = 0.5; 
figure(2); hold on;
for i = 1:m
    nu = nu/2; K = (nu/h/h)*A + (b/h)*B;  
    u = K\ones(N,1); uu = [0;u;0];
    plot(xx,uu,'-','LineWidth',1.5);
end
xlabel('x'); ylabel('u'); grid on; saveas(2,'cd1-uwd.pdf');
end
%%%%