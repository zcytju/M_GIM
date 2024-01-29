%----------------DCB and IONC estimate based on multi stations----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear all;
%---------------------------------Setting--------------------------------------
%--Rinex files
r_ipath='Files_rinex\global';
%--SP3 files
s_ipath='Files_SP3';
%--DCB files
i_ipath='Files_dcb';
%--system
% G : GPS-only    GC : G+C    GE : G+E    GR : G+R    GCE: G+C+E
% GCR: G+C+R    GER: G+E+R    GCER : G+C+E+R
system='GCER';
%--cut-off angle
lim=10;
%--number of ionosphere model parameter groups
fig=12;
% --sampling interval
sample_interval=30;
%--SH model order and degree
order=15;
%--the weight of each system observation
PG=1; PC=1/2; PE=1; PR=1/3;
%--the cores number of computer processor
Corenum=2;

tic
%Step one--------------------Read RINEX files------------------------------
disp('M_GIM starts running!')
global sample_num;
sample_num=24*60^2/sample_interval;
% --L1,L2,P1,P2 and coordinates of GPS receviers are obtained by this step
r_opath='OBS\global';
disp('Step one: read rinex files !')
Sites_Info=read_rinex(r_ipath,r_opath);
disp('Step one: completing !')
%Step two--------------------Read SP3 files--------------------------------
%--Get satellites coordinates and interpolation
%--SP3 files of the previous day and the next day are needed for interpolation
s_opath='SP3';
disp('Step two: read SP3 files !')
load('sate_mark.mat');
read_sp3(s_ipath,s_opath,sate_mark);
disp('Step two: completing !')
%Step three-------------------Read DCB files------------------------------
%--Get reference information
disp('Step three: read DCB files !')
[SDCB_REF,Sites_Info]= read_dcb(i_ipath,Sites_Info,sate_mark);
disp('Step three: completing !')
%
%Step four---------------Ionosphere Observations-----------------------
%--get smoothed P4 observations
disp('Step four: data pre-processing and get ionosphere observations !')
p_opath='P4\global';
addpath('Func_SH','Tools\M_DCB');
Get_P4(r_opath,p_opath,s_opath,Sites_Info,lim*pi/180,sate_mark);
disp('Step four: completing !')
%Step five-----------------------Get M_GIM---------------------------------
list_sp3=dir('SP3\*.mat');
list=dir('P4\global\GPS');
len=length(list);
disp('Step five: IONC and DCB estimate day by day!')
result_folder='M_Result';
if ~exist(result_folder,'dir')
    mkdir(result_folder);
end
if license('test', 'Distrib_Computing_Toolbox')
    mypool=parpool('local',Corenum);
end
%--chose the order of spheric harmonic function
for i=3:len
    doy=list(i).name;
    total_station=length(dir(['OBS\global\',doy,'\*.mat']));
    index=find_sp3(list_sp3,doy);
    load(['SP3\' list_sp3(index).name],'-mat');
    list_P4.gps=dir([p_opath '\GPS\',doy,'\*.mat']);
    list_P4.bds=dir([p_opath '\BDS\',doy,'\*.mat']);
    list_P4.gal=dir([p_opath '\GAL\',doy,'\*.mat']);
    list_P4.galx=dir([p_opath '\GALX\',doy,'\*.mat']);
    list_P4.glo=dir([p_opath '\GLO\',doy,'\*.mat']);
    
    if strcmp(system,'G')
        [G_R, G_S,IONC,m0,NN]= Get_SH_G(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,sate_mark);
        fname=['M_Result\G' doy '.mat'];
        save(fname,'G_R','G_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GC')
        [G_R, G_S, C_R, C_S, IONC,m0,NN]= Get_SH_GC(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PC,sate_mark);
        fname=['M_Result\GC' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GE')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, IONC,m0,NN]= Get_SH_GE(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PE,sate_mark);
        fname=['M_Result\GE' doy '.mat'];
        save(fname,'G_R','G_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GR')
        [G_R, G_S, R_R, R_S, IONC,m0,NN]= Get_SH_GR(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PR,sate_mark);
        fname=['M_Result\GR' doy '.mat'];
        save(fname,'G_R','G_S','R_R','R_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GCE')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, C_R, C_S, IONC,m0,NN]= Get_SH_GCE(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PC,PE,sate_mark);
        fname=['M_Result\GCE' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GCR')
        [G_R, G_S, C_R, C_S, R_R, R_S, IONC,m0,NN]= Get_SH_GCR(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PC,PR,sate_mark);
        fname=['M_Result\GCR' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','R_R','R_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GER')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, R_R, R_S, IONC,m0,NN]= Get_SH_GER(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PE,PR,sate_mark);
        fname=['M_Result\GER' doy '.mat'];
        save(fname,'G_R','G_S','R_R','R_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GCER')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, C_R, C_S, R_R, R_S, IONC,m0,NN]= Get_SH_GCER(fig,doy,Sites_Info,sate,SDCB_REF,order,PG,PC,PE,PR,sate_mark);
        fname=['M_Result\GCER' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','R_R','R_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    clear sate;
    disp(['DCB and IONC estimate of doy ' doy ' complete!']);
end
if license('test', 'Distrib_Computing_Toolbox')
    delete(mypool);
end
disp('Step five: completing!')
toc