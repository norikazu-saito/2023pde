function plotTriangulation(datafile)
% get mesh data
[deg,p,idp,t,idt,e,ide,d,idd] = getmesh(datafile);
% triangulation 
figure(2); 
%triplot(t',p(1,:),p(2,:),'cyan','LineWidth',2); 
triplot(t',p(1,:),p(2,:),'blue','LineWidth',1); 
%
ax = gca; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; 
axis equal; saveas(2,"triangulation1.pdf");
end