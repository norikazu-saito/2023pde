function nschr(N, theta, Tmax)
    % parameters
    a = 0; b = 2*pi; h = (b-a)/(N+1); x = (a+h:h:b-h).'; xx = (a:h:b).';
    ua = 0 + i*0; ub = 0 + i*0;
    tau = 0.1*h; lambda = tau/h/h; nmax = floor(Tmax/tau);
    number = 40; step = floor(max(1,nmax/number));
    % discrete Laplacian
    MM = 2*eye(N,N) - diag(ones(N-1,1),-1) - diag(ones(N-1,1),1);
    H = sparse(eye(N,N) + 1i * theta * lambda * MM);
    K = sparse(eye(N,N) - 1i * (1-theta) * lambda * MM);
    [L,U,P] = lu(H);
    % t = 0
    tnow = 0; u = func_a(x); uu = [ua;u;ub];
    % for plotting
    figure(1); hold on; tsp = tnow*ones(1,N+2); plot3(xx,tsp,abs(uu),'r');
    % t > 0
    kmax = 2000; tol = 1e-8;
    for n = 1:nmax
        tpast = tnow; tnow = n*tau; g0 = func_g(tpast,x); g1 = func_g(tnow,x);
        % relaxation iteration
        k = 1; dif = tol + 1; w = u;
        while (k < kmax) && (dif > tol)
            f0 = (1-theta)*(func_f(u)+g0); f1 = theta*(func_f(w)+g1);
            f = H*w - K*u + 1i*tau*(f0+f1);
            v = (U\(L\(P*f))); w = w - v; dif = norm(v,inf); k = k + 1;
        end
        if k == kmax
            disp("relaxation iteration does not converge!");
            disp(tnow); break;
        end
        u = w;
        if rem(n,step) == 0
            % for plotting
            figure(1); tsp = tnow*ones(1,N+2); uu = [ua;u;ub];
            plot3(xx,tsp,abs(uu),'b');
        end
    end
    % decoration
    figure(1); xlabel('x'); ylabel('t'); zlabel('abs(u)'); grid on;
    view(60,52); saveas(1,'nschr.pdf')
end
% initial value
function w = func_a(x)
    w = sin(x).^2 + i*sin(x).^2;
end
% inhomogeneous term
function w = func_g(t,x)
    w = 0 + i*0;
end
% nonlinearity
function w = func_f(z)
    % w = z.*abs(z).^2;
    w = z.*abs(z).^(0.5);
end
