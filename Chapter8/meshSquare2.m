function meshSquare2(L, m)
h=L/m; mm=m+1; Nn=mm^2; Ne=2*m^2; Nr=0; Nd=4*m; 
x=linspace(0,L,mm); dat=[]; t=[]; p=[]; er=[]; ed=[]; 
% total numbers 
dat=[dat;Nn,Ne,Nr,Nd];
% nodes
for j=1:mm
    for i=1:mm
        k=(j-1)*mm + i; id=0; 
        tmp=(i-1)*(j-1)*(i-mm)*(j-mm);
        if tmp==0 id=1; end
        p=[p;x(i),x(j),id,0]; 
    end
end
% elements
enum=0; 
for i=1:m
 for j=1:m
   k=(j-1)*mm + i; kk=k+mm;
   % below
   t = [t;k, k+1, kk+1, 0]; 
   % above
   t=[t;k, kk+1, kk, 0]; 
 end
end
% Dirichlet boundaries
for j=1:m
ed=[ed;j,j+1,1,0];
end
for j=1:m
 q=j*mm; ed=[ed;q,q+mm,1,0];
end
for i=1:m
 q = mm^2-i+1; ed=[ed;q,q-1,1,0];
end
for j=1:m
 q=mm^2-m-(j-1)*mm;ed=[ed;q,q-mm,1,0];
end
% output 
dat=[dat;p;t;er;ed]';
F1 = fopen("square2.dat","w"); fprintf(F1,"%f %f %f %f\n",dat); fclose(F1);
end