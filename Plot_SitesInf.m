%% plot station distributions with quad-system results
%--19: year; 275: day
figure;    doy=19275; 
% Set latitude and longitude range
lat2=-90;       lat1=90;
lon1=-180;    lon2=180;
load(['M_Result/GCER',num2str(doy),'.mat']);
warning off;
addpath('Tools/m_map','Tools/m_map/private');
[site_inf]=Plot_Multisites(lat2,lat1,lon1,lon2,G_R,C_R,E_R,EX_R,R_R,Sites_Info,list_P4);
saveas(gcf,['M_PLOT/','Global stations.fig'],'fig');

%% --------------------------------------------------------------------------------
function [Site_inf]=Plot_Multisites(lat1,lat2,lon1,lon2,G_R,C_R,E_R,EX_R,R_R,Sites_Info,list_P4)
%% plot station distributions with quad-system results
% INPUT:
%     lat1, lat2, lon1, lon2: latitude and longitude range
%     G_R, C_R, E_R, EX_R, R_R: receiver DCBs of different systems
%     Sites_Info: name and coordinate information of the stations
%     list_P4: P4 observation list for different systems
% OUTPUT:
%     Site_inf: reconstructed station information
%% written by Zhou C. et al., 2020/7/15
%% --------------------------------------------------------------------------------
list_G=list_P4.gps;
G_r=length(list_G);%the number of receivers
tol_G=size(Sites_Info.name,2);
inf_G=[];
for i=1:G_r
    site=list_G(i).name(1:4);
    G_real=G_R(i);
    for j=1:tol_G
        if string(site)==string(Sites_Info.name(j))
            [latitude, longitude]=XYZtoBLH(Sites_Info.coor(j,1),Sites_Info.coor(j,2),Sites_Info.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            RDCB_ref=Sites_Info.GPS_RDCB_REF(j);
        end
    end
    inf_G=[inf_G;string(site), latitude, longitude, G_real, RDCB_ref];
end
Site_inf.inf_G=inf_G;
disp('------> GPS data has been constructed!');

list_R=list_P4.glo;
R_r=length(list_R);%the number of receivers
tol_R=size(Sites_Info.name,2);
inf_R=[];
for i=1:R_r
    site=list_R(i).name(1:4);
    R_real=R_R(i);
    for j=1:tol_R
        if string(site)==string(Sites_Info.name(j))
            [latitude, longitude]=XYZtoBLH(Sites_Info.coor(j,1),Sites_Info.coor(j,2),Sites_Info.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            RDCB_ref=Sites_Info.GLO_RDCB_REF(j);
        end
    end
    inf_R=[inf_R;string(site), latitude, longitude, R_real, RDCB_ref];
end
Site_inf.inf_R=inf_R;
disp('------> GLONASS data has been constructed!');

list_E=list_P4.gal;
E_r=length(list_E);%the number of receivers
tol_E=size(Sites_Info.name,2);
inf_E=[];
for i=1:E_r
    site=list_E(i).name(1:4);
    E_real=E_R(i);
    for j=1:tol_E
        if string(site)==string(Sites_Info.name(j))
            [latitude, longitude]=XYZtoBLH(Sites_Info.coor(j,1),Sites_Info.coor(j,2),Sites_Info.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            RDCB_ref=Sites_Info.GAL_RDCB_REF(j);
        end
    end
    inf_E=[inf_E;string(site), latitude, longitude, E_real, RDCB_ref];
end
Site_inf.inf_E=inf_E;

list_EX=list_P4.galx;
EX_r=length(list_EX);%the number of receivers
tol_EX=size(Sites_Info.name,2);
inf_EX=[];
for i=1:EX_r
    site=list_EX(i).name(1:4);
    EX_real=EX_R(i);
    for j=1:tol_EX
        if string(site)==string(Sites_Info.name(j))
            [latitude, longitude]=XYZtoBLH(Sites_Info.coor(j,1),Sites_Info.coor(j,2),Sites_Info.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            RDCB_ref=Sites_Info.GALX_RDCB_REF(j);
        end
    end
    inf_EX=[inf_EX;string(site), latitude, longitude, G_real, RDCB_ref];
end
Site_inf.inf_EX=inf_EX;
disp('------> Galileo data has been constructed!');

list_C=list_P4.bds;
C_r=length(list_C);%the number of receivers
tol_C=size(Sites_Info.name,2);
inf_C=[];
for i=1:C_r
    site=list_C(i).name(1:4);
    C_real=C_R(i);
    for j=1:tol_C
        if string(site)==string(Sites_Info.name(j))
            [latitude, longitude]=XYZtoBLH(Sites_Info.coor(j,1),Sites_Info.coor(j,2),Sites_Info.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            RDCB_ref=Sites_Info.BDS_RDCB_REF(j);
        end
    end
    inf_C=[inf_C;string(site), latitude, longitude, C_real, RDCB_ref];
end
Site_inf.inf_C=inf_C;
disp('------> BDS data has been constructed!');

scatter(1, 1,'Marker','o','MarkerEdgeColor',[0.00,0.00,1.00],...
        'MarkerFaceColor',[0.00,0.00,1.00], 'LineWidth', 5);  % plot the gps sites
    hold on;
scatter(1, 2,'Marker','square','MarkerEdgeColor',[0.93,0.69,0.13],...
        'MarkerFaceColor',[0.93,0.69,0.13], 'LineWidth', 3.5);  % plot the glonass sites
    hold on;
scatter(1, 3,'Marker','o','MarkerEdgeColor',[0.47,0.67,0.19],...
        'MarkerFaceColor',[0.47,0.67,0.19], 'LineWidth', 1.8);  % plot the galileo sites
    hold on;
scatter(1, 4,'Marker','square','MarkerEdgeColor',[0.64,0.08,0.18],...
        'MarkerFaceColor',[0.64,0.08,0.18], 'LineWidth', 0.2);  % plot the bds sites
    hold on;   

% m_proj('Miller Cylindrical','longitudes',[lon1 lon2], 'latitudes',[lat1,lat2]);
m_proj('Robinson','longitudes',[lon1 lon2], 'latitudes',[lat1,lat2]);
m_coast('patch',[0.92,0.92,0.92]);
m_grid('xaxis','bottom','FontSize',12);%'box','fancy',    ,'FontSize',14,'FontName','Times New Roman'
hold on;

G_num=size(inf_G,1);
for i=1:G_num
    lat=str2double(inf_G(i,2));
    lon=str2double(inf_G(i,3));
    m_scatter(lon, lat,'Marker','o','MarkerEdgeColor',[0.00,0.00,1.00],...
        'MarkerFaceColor',[0.00,0.00,1.00], 'LineWidth', 5)  % plot the gps sites
    hold on;
%     m_text(lon, lat, inf_G(i,1),'fontsize',12);  % marking the sites,'FontName','Times New Roman'
%     hold on;
end
disp('------> GPS data has been plotted!');

R_num=size(inf_R,1);
for i=1:R_num
    lat=str2double(inf_R(i,2));
    lon=str2double(inf_R(i,3));
    m_scatter(lon, lat,'Marker','square','MarkerEdgeColor',[0.93,0.69,0.13],...
        'MarkerFaceColor',[0.93,0.69,0.13], 'LineWidth', 3.5)  % plot the sites
    hold on;
end
disp('------> GLONASS data has been plotted!');

E_num=size(inf_E,1);
for i=1:E_num
    lat=str2double(inf_E(i,2));
    lon=str2double(inf_E(i,3));
    m_scatter(lon, lat,'Marker','o','MarkerEdgeColor',[0.47,0.67,0.19],...
        'MarkerFaceColor',[0.47,0.67,0.19], 'LineWidth', 1.8)  % plot the sites
    hold on;
end

EX_num=size(inf_EX,1);
for i=1:EX_num
    lat=str2double(inf_EX(i,2));
    lon=str2double(inf_EX(i,3));
    m_scatter(lon, lat,'Marker','o','MarkerEdgeColor',[0.47,0.67,0.19],...
        'MarkerFaceColor',[0.47,0.67,0.19], 'LineWidth', 1.8)  % plot the sites
    hold on;
end
disp('------> Galileo data has been plotted!');

C_num=size(inf_C,1);
for i=1:C_num
    lat=str2double(inf_C(i,2));
    lon=str2double(inf_C(i,3));
    m_scatter(lon, lat,'Marker','square','MarkerEdgeColor',[0.64,0.08,0.18],...
        'MarkerFaceColor',[0.64,0.08,0.18], 'LineWidth', 0.2)  % plot the sites
    hold on;
end
disp('------> BDS data has been plotted!');
legend('','GPS','GLONASS','Galileo','BDS','Location','northwest','Orientation','horizontal' ,'FontSize',14,'Box','off');    
end

