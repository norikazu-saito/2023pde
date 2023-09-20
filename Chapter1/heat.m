function [h, error1, error2] = heat(N, theta, lam, Tmax, givenfunc, id)
    % id >= 0 then plot graph;  id < 0  then no graph
    % parameters
    k = 1; a = 0; b = 1; h = (b-a)/(N+1); x = (a+h:h:b-h).'; xx = (a:h:b).';
    ua = 0; ub = 0; tau = lam*h^2/k; nmax = floor(Tmax/tau);
    number = 40; step = floor(max(1,nmax/number));
    % discrete Laplacian
    A = 2*eye(N,N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
    H = sparse(eye(N,N) + theta*lam*A); 
    K = sparse(eye(N,N) - (1-theta)*lam*A); 
    [L,U,P] = lu(H);
    % t = 0
    tnow = 0.0; u = givenfunc(0,x,tnow); uu = [ua;u;ub];
    % for plotting
    if id >= 0
        figure(1); hold on; plot(xx,uu,'r');
        figure(2); hold on; tsp = tnow*ones(1,N+2); plot3(xx,tsp,uu,'r');
    end
    % t > 0
    error1 = -1; error2 = -1;
    for n = 1:nmax
        tpast = tnow; tnow = n*tau;
        f = K*u + tau*((1-theta)*givenfunc(1,x,tpast) + ...
            theta*givenfunc(1,x,tnow));
        u = (U\(L\(P*f))); uu = [ua;u;ub];
        % error
        err = u - givenfunc(2,x,tnow); error1 = max(error1,norm(err,inf));
        error2 = max(error2,norm(err,2)); error2 = error2*sqrt(h);
        % for plotting
        if id >= 0
            if rem(n,step) == 0
                figure(1); plot(xx,uu,'b');
                figure(2); tsp = tnow*ones(1,N+2); plot3(xx,tsp,uu,'b');
            end
        end
    end
    if id >= 0
        % decoration
        figure(1); xlabel('x'); ylabel('u'); grid on; saveas(1,'heatA.pdf');
        figure(2); xlabel('x'); ylabel('t'); zlabel('u'); grid on;
        view(60,15); saveas(2,'heatB.pdf');
    end
end
