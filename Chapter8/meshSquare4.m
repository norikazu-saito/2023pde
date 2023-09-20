function meshSquare4(L, m)
h=L/m; mm=m+1; np=mm^2; ne=2*m^2; nr=4*m; nd=0;
x=linspace(0,L,mm); dat=[]; t=[]; p=[]; er=[]; ed=[]; 
% total numbers 
dat=[dat;np,ne,nr,nd];
% nodes
for j=1:mm
    for i=1:mm
        k=(j-1)*mm + i; id=0; tmp=(i-1)*(j-1)*(i-mm)*(j-mm);
        if tmp==0 id=1; end 
        p=[p;x(i),x(j),id,0]; 
    end
end
% elements
enum=0; 
for j=1:m
    for i=1:m
        k=(j-1)*mm + i; kk=k+mm;
        t = [t;k, k+1, kk+1, 0]; t=[t;k, kk+1, kk, 0]; 
    end
end
% neumann boundaries
for i=1:m er=[er;i,i+1,2,0]; end
for j=1:m q=j*mm; ed=[ed;q,q+mm,2,0]; end
for i=1:m q=mm^2-i+1; ed=[ed;q,q-1,2,0]; end
for j=1:m q=mm^2-m-(j-1)*mm;ed=[ed;q,q-mm,2,0]; end
% output 
dat=[dat;p;t;er;ed]'; F1 = fopen("square4.dat","w"); 
fprintf(F1,"%f %f %f %f\n",dat); fclose(F1);
end