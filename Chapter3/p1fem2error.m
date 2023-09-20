function p1fem2error
% preliminaries 
hvec=[]; l2vec=[]; h1vec=[]; linfvec=[]; hv1=[]; hv2=[]; 
%str = ["square2010.dat","square2020.dat","square2040.dat","square2080.dat","square2160.dat"];
str = ["square2010.dat","square2020.dat","square2040.dat"];
% 
for it=1:size(str,2)
    [hsize,l2,h1,linf] = p1fem2(str(it),@femFunc2,1);
    hvec=[hvec;hsize]; l2vec=[l2vec;l2]; 
    h1vec=[h1vec;h1]; linfvec=[linfvec;linf];
    if it==1
        sl1 = 1.5*max(h1,l2); sl2 = 0.5*min(h1,l2); 
        hv1 = [hv1;sl1]; hv2 = [hv2;sl2]; 
    else 
        sl1=sl1/2; sl2=sl2/4; hv1 = [hv1;sl1]; hv2 = [hv2;sl2];
    end
end
%
hvec1 = log(hvec)-log(circshift(hvec,1)); 
l2vec1 = log(l2vec)-log(circshift(l2vec,1));
h1vec1 = log(h1vec)-log(circshift(h1vec,1)); 
linfvec1 = log(linfvec)-log(circshift(linfvec,1)); 
ratel2=l2vec1./hvec1; ratel2(1)=0; rateh1=h1vec1./hvec1; rateh1(1)=0;
ratelinf = linfvec1 ./ hvec1; ratelinf(1)=0;
Z=[hvec';l2vec';h1vec';linfvec';ratel2';rateh1';ratelinf'];
% output resuts
F1=fopen("p1fem2error.res","w"); 
fprintf(F1,"%f %.5e %.5e %.5e %.4f %.4f %.4f\n",Z); fclose(F1);
% plot
figure(7); 
loglog(hvec,l2vec,'-xr','LineWidth',2,'MarkerSize',10); hold on;
loglog(hvec,linfvec,':ob','LineWidth',2,'MarkerSize',10);
loglog(hvec,h1vec,'-^m','LineWidth',2,'MarkerSize',10); hold off; 
% decoration
legend('L^2 err','L^{inf} err','H^1 err','Location','SouthEast');
xlabel('h');ylabel('ERROR'); grid on; pbaspect([1 2 1]); 
xticks([1e-2 1e-1]); saveas(7,'p1fem2error.pdf'); 
end