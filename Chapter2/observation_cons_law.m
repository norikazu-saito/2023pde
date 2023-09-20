function observation_cons_law(N, Tmax)
[res05]=lschr11(N, 0.5, Tmax); [res10]=lschr11(N, 1, Tmax);
time=res05(:,1); mass05=res05(:,2); ener05=res05(:,3);
mass10=res10(:,2); ener10=res10(:,3);
%
figure(5); hold on; 
plot(time,mass05,'-r','LineWidth',4); 
plot(time,mass10,':b','LineWidth',4); hold off; 
figure(7); hold on; 
plot(time,ener05,'-r','LineWidth',4); 
plot(time,ener10,':b','LineWidth',4); hold off; 
%
figure(5);xlabel('time');ylabel('mass'); grid on; 
legend('theta: 0.5','theta: 1','Location','SouthWest');
figure(7);xlabel('time');ylabel('energy'); grid on; 
legend('theta: 0.5','theta: 1','Location','SouthWest');
saveas(5,'conslaw05.pdf'); saveas(7,'conslaw10.pdf');   
end