function p1femBT(datafile,givenfunc,Tmax)
% read mesh data
[deg,p,idp,t,idt,e,ide,d,idd] = getmesh(datafile);
% diffusion coef
coef_d=0.01;
% matrices
[K, M, C, kappa, bmax] = matrixBT(p,t,givenfunc); 
Id = eye(size(K)); B = diag(1./diag(M))*(coef_d*K + C);
% initial value
usol=givenfunc(1,p(1,:),p(2,:))'; 
fign=1; figure(fign); plotP1fem(p,t,usol); 
umass=[]; umass=[umass;dot(M*usol,ones(size(usol)))]; 
time=[]; time=[time;0];
%% time increment  
tau = 0.9*kappa^2/(3*coef_d + 4*kappa*bmax); nmax=floor(Tmax/tau); 
number = 6; step = floor(max(1,nmax/number)); ptime=[];ptime=[ptime;0];
for n=1:nmax   
    tnow = n*tau; 
    usol = (Id - tau*B) * usol;
    umass=[umass;dot(M*usol,ones(size(usol)))]; time=[time;tnow]; 
    %% for 2d & 3d draw
    if rem(n, step)==0
        fign=fign+1; figure(fign); plotP1fem(p,t,usol); 
        ptime=[ptime;tnow];
        fname = "p1fem1convect" + fign + ".pdf"; saveas(fign,fname);
    end
end
% show results
for n=1:fign
    fname = "p1femBT" + n + ".pdf"; 
    figure(n);xlabel('x');ylabel('y');zlabel('u');grid on;
    title("time="+ptime(n)); 
    view(45,45);saveas(n,fname);
end
% mass conservation
fign=fign+10; figure(fign); plot(time,umass,'LineWidth',2);
xlabel('time');ylabel('mass');saveas(fign,"p1BTm.pdf");
% important values
kappa
tau
%%% end of program
end