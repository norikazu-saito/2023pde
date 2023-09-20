function plotP1fem(p, t, u)
    nt = size(t,2); X = zeros(3*nt,1); Y = zeros(3*nt,1); uu = zeros(3*nt,1);
    i = t(1,:); j = t(2,:); k = t(3,:);
    X(1:3:end) = p(1,i); X(2:3:end) = p(1,j); X(3:3:end) = p(1,k);
    Y(1:3:end) = p(2,i); Y(2:3:end) = p(2,j); Y(3:3:end) = p(2,k);
    uu(1:3:end) = u(i); uu(2:3:end) = u(j); uu(3:3:end) = u(k);
    % using full color
    trisurf(reshape([1:3*nt],3,nt)',X,Y,uu);
    view(45,30);
    % axes settings
    xlabel('x','FontSize',14); ylabel('y','FontSize',14); 
end
