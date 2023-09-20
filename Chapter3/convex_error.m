function convex_error
formatSpec = '%f %f %f %f %f';sizeD = [5 Inf]; 
FILE1 = fopen("nonconvex1.dat",'r'); 
FILE2 = fopen("nonconvex2.dat",'r'); 
FILE3 = fopen("nonconvex3.dat",'r'); 
D1 = fscanf(FILE1,formatSpec,sizeD); 
D2 = fscanf(FILE2,formatSpec,sizeD); 
D3 = fscanf(FILE3,formatSpec,sizeD); 
fclose(FILE1); fclose(FILE2); fclose(FILE3);
% D1(1,:)でsize D1(2,:)でH1error
hvec=D1(1,:);vec1=D1(2,:);vec2=D2(2,:);vec3=D3(2,:);
% 
hvec5=[D1(1,3),D1(1,6)]; f1=0.03;
vec5=[f1,exp(log(f1)-1*(log(hvec5(1)/hvec5(2))))];
hvec6=[D1(1,3),D1(1,6)]; f1=0.2;
vec6=[f1,exp(log(f1)-0.5*(log(hvec5(1)/hvec5(2))))];
% plot
figure(6); 
loglog(hvec,vec1,'-sr','LineWidth',2,'MarkerSize',10); hold on;
loglog(hvec,vec2,'-xb','LineWidth',2,'MarkerSize',10); 
loglog(hvec,vec3,':om','LineWidth',2,'MarkerSize',10);
loglog(hvec5,vec5,'-k','LineWidth',2,'MarkerSize',10);
loglog(hvec6,vec6,'--k','LineWidth',2,'MarkerSize',10);
hold off;
% decoration
legend('0.7pi','1.5pi','1.8pi','slope 1','slope 1/2','Location','NorthWest');
xlabel('h');ylabel('ERROR'); grid on; pbaspect([1 2 1]); 
xticks([1e-2 1e-1]); saveas(6,'nonconvec10error.pdf'); 
end