%%% initial value
function w = init1(x)
    w = x .* (sin(3*pi*x)).^2; 
end