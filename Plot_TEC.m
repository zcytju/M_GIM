function Plot_TEC(fig,latlim,lonlim,Pname,PaintData,lat1,lat2,lon1,lon2,limmin,limmax)
%% plot global or regional ionospheric maps
% INPUT:
%     fig: group number of products
%     latlim, lonlim: latitude and longitude resolution
%     Pname: specified name
%     PaintData: data need to be plotted
%     lat2, lat1: latitude range
%     lon1, lon2: longitude range
%     limmin, limmax: lower and upper limits of the bar value
%% written by Zhou C. et al., 2020/7/15
%% ============================================
warning off;
result_folder='M_PLOT';
if ~exist(result_folder,'dir')
    mkdir(result_folder);
end
num=((lat1-lat2)/latlim+1)*((lon2-lon1)/lonlim+1);
for i=1:fig
	x=PaintData(1+num*(i-1):num*i,3); 
	y=PaintData(1+num*(i-1):num*i,2); 
	z=PaintData(1+num*(i-1):num*i,4); 
	[X,Y,Z]=griddata(x,y,z,linspace(lon1,lon2,101)',linspace(lat2,lat1,101),'v4'); 
    contourf(X,Y,Z);
    lat=lat2:((lat1-lat2)/100):lat1;
    lon=lon1:((lon2-lon1)/100):lon2;
    [Plg,Plt]=meshgrid(lon,lat); %constract grid
    m_proj('Miller Cylindrical','longitudes',[lon1 lon2], 'latitudes',[lat2 lat1]);
    m_pcolor(Plg,Plt,Z);
    set(gca,'FontSize',14);%,'FontName','Times New Roman'
    hold on;
    m_coast('line','color','k','linewidth',0.5);
    m_grid('xaxis','bottom');
	shading interp;
	colormap(jet);
	h=colorbar;
    set(get(h,'Title'),'string','TECU');
	caxis([limmin,limmax]);
	flag(i)=(i-1)*(24/fig);
    set(gca,'FontSize',14);%,'FontName','Times New Roman'
	title([Pname,num2str(flag(i)),' : 00'],'FontSize',16);
	axis tight;
	set(gca,'nextplot','replacechildren');
	F(i)=getframe;
	picname=[Pname,num2str(i),'.fig'];
 	saveas(gcf,['M_PLOT\',picname]);
    disp(['--------> [ ',num2str(i),' / ',num2str(fig),' ] ionospheric maps have been plotted!']);
end
close all;
% movie(F,fig);
%
% for i=1:fig
%     picname=[Pname,num2str(i),'.fig'];
%     open(['M_PLOT\',picname]);
%     frame=getframe(gcf);  
%     im=frame2im(frame); 
%     [I,map]=rgb2ind(im,20); 
%     filename=['M_PLOT\',Pname,'.gif'];
%     if i==1
%         imwrite(I,map,filename,'gif', 'Loopcount',inf,'DelayTime',0.5);
%     elseif i==fig
%         imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.5);
%     else
%         imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.5);
%     end  
%     close all;
% end
% disp('--------> Ionospheric maps in .gif format have been plotted!');

end
