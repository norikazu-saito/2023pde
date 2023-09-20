function [w, Q] = gauss2D(x, y, area)
    % weights
    d = [(155-sqrt(15))/1200; (155+sqrt(15))/1200];
    w = [9/40;d(1);d(1);d(1);d(2);d(2);d(2)]; w = w*area;
    % quadratue points
    v = [(6-sqrt(15))/21; (6+sqrt(15))/21];
    theta = [1/3,1/3,1/3;
             v(1),v(1),1-2*v(1); v(1),1-2*v(1),v(1); 1-2*v(1),v(1),v(1);
             v(2),v(2),1-2*v(2); v(2),1-2*v(2),v(2); 1-2*v(2),v(2),v(2)];
    Z = [x',y']; Q = theta*Z;
end
