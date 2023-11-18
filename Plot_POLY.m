%% ==========Plot and write regional ionospheric maps=================
%% polynomial model
doy=19275;  fig=24;  K=6;  M=4;
load('sate_mark.mat');
load(['M_Result/POLY_GCER',num2str(doy),'.mat']);
warning off;
addpath('Tools/m_map','Tools/m_map/private');
lat2=40;     lat1=80;    lon1=0;    lon2=60;
latlim=2.5;   lonlim=5;
lat0=deg2rad(60);  lon0=deg2rad(30);
VVTEC=Get_VTEC_POLY(fig,latlim,lonlim,IONC,m0,NN,K,M,lat2,lat1,lon1,lon2,lat0,lon0);
VTEC=VVTEC;VTEC(VTEC(:,4)<0,4)=0.05;
RMS=[VTEC(:,1:3),VTEC(:,5)];
% % read CODE final GIMs (codg2750.19i) as reference
disp('--------> Read CODE final GIMs as reference !');
IGSData=read_ionex(fig,'TEC');
AreaTEC=Get_areaTEC(fig,lat2,lat1,lon1,lon2,IGSData);
DIFFTEC=[VTEC(:,1:3),VTEC(:,4)-AreaTEC(1:size(VTEC,1),4)];
Pname='2019-275-POLY- '; 
Plot_TEC(fig,latlim,lonlim,Pname,VTEC,lat1,lat2,lon1,lon2,0,15);
Pname_D='2019-275-POLY-D- '; 
Plot_TEC(fig,latlim,lonlim,Pname_D,DIFFTEC,lat1,lat2,lon1,lon2,-5,5);
Pname_R='2019-275-POLY-RMS- '; 
Plot_TEC(fig,latlim,lonlim,Pname_R,RMS,lat1,lat2,lon1,lon2,0,1);
% write results to ionex files
Write_ionex(doy,fig,G_R,C_R,E_R,EX_R,R_R,G_S,C_S,E_S,R_S,Pname,VTEC,sate_mark,total_station,list_P4);
%% ++++++++++++++++PLOT OVER!!!+++++++++++++++++++