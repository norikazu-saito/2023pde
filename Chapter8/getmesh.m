function [deg, p, idp, t, idt, e, ide, d, idd] = getmesh(filename)
    FILE1 = fopen(filename,'r'); formatSpec = '%f %f %f %f';
    sizeD = [4 Inf]; D = fscanf(FILE1, formatSpec, sizeD); fclose(FILE1);
    D = D'; np = D(1,1); nt = D(1,2); ne = D(1,3); nd = D(1,4);
    deg = [np nt ne nd]; p = D(2:np+1, 1:2)'; idp = D(2:np+1, 3)';
    t = D(np+2:np+nt+1, 1:3)'; idt = D(np+2:np+nt+1, 4)';
    e = D(np+nt+2:np+nt+ne+1, 1:2)'; ide = D(np+nt+2:np+nt+ne+1, 3)';
    d = D(np+nt+ne+2:np+nt+ne+nd+1, 1:2)';
    idd = D(np+nt+ne+2:np+nt+ne+nd+1, 3)';
end
