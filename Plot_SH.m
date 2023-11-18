%% ----------------- Plot and write global ionospheric VTEC maps -------------------
%% =================================================
doy=19275;  fig=12; 
load('sate_mark.mat');
load(['M_Result/GCER',num2str(doy),'.mat']);
warning off;
addpath('Tools/m_map','Tools/m_map/private');
lat2=-87.5;    lat1=87.5; lon1=-180;    lon2=180;
VVTEC = Get_VTEC(fig, 2.5, 5, IONC, NN, m0, 15);
VTEC=VVTEC;VTEC(VTEC(:,4)<0,4)=0.05;
RMS=[VTEC(:,1:3),VTEC(:,5)];
% % read IGS final GIMs (igsg2750.19i) as reference
disp('--------> Read IGS final GIMs as reference !');
IGSData=read_ionex(fig,'TEC');
DIFFTEC=[VTEC(:,1:3),VTEC(:,4)-IGSData(1:size(VTEC,1),4)];
Pname='2019-275-GCER- '; 
Plot_TEC(fig,2.5,5,Pname,VTEC,lat1,lat2,lon1,lon2,0,50);
Pname_D='2019-275-GCER-D- '; 
Plot_TEC(fig,2.5,5,Pname_D,DIFFTEC,lat1,lat2,lon1,lon2,-10,10);
Pname_R='2019-275-GCER-RMS- '; 
Plot_TEC(fig,2.5,5,Pname_R,RMS,lat1,lat2,lon1,lon2,0,5);
% write results to ionex files
Write_ionex(doy,fig,G_R,C_R,E_R,EX_R,R_R,G_S,C_S,E_S,R_S,Pname,VTEC,sate_mark,total_station,list_P4);

%% ++++++++++++++++PLOT OVER!!!+++++++++++++++++++