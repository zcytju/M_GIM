%----------------DCB and IONC estimate based on multi stations----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear all;
%---------------------------------Setting-------------------------------------- 
%--Rinex files
r_ipath='Files_rinex\regional';
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
fig=24;
%--nonintergral SH model order and degree
K=6; M=4;
%--the weight of each system observation
PG=1; PC=1/2; PE=1; PR=1/3;
%--the cores number of computer processor
Corenum=2;

tic
%Step one--------------------Read RINEX files------------------------------
disp('MDCB(multi-stations) starts running!');
% --L1,L2,P1,P2 and coordinates of GPS receviers are obtained by this step
r_opath='OBS\regional';
disp('Step one: read rinex files !');
% Sites_Info=read_rinex(r_ipath,r_opath);
disp('Step one: completing !');
%Step two--------------------Read SP3 files--------------------------------
%--Get satellites coordinates and interpolation
%--SP3 files of the previous day and the next day are needed for interpolation
s_opath='SP3';
disp('Step two: read SP3 files !');
load('sate_mark.mat');
% read_sp3(s_ipath,s_opath,sate_mark);
disp('Step two: completing !');
%Step three------------------Read ionex files------------------------------
%--Get reference information
disp('Step three: read DCB files !');
% [SDCB_REF,Sites_Info]= read_dcb(i_ipath,Sites_Info,sate_mark);
disp('Step three: completing !');
% 
%Step four---------------Ionosphere Observations-----------------------
%--get smoothed P4 observations
disp('Step four: data pre-processing and get ionospheric observations !')
p_opath='P4\regional';
addpath('Func_nonSH','Tools\M_DCB');
% Get_P4(r_opath,p_opath,s_opath,Sites_Info,lim*pi/180,sate_mark);
disp('Step four: completing !')
%Step five----------------------Get M_GIM----------------------------------
list_sp3=dir('SP3\*.mat');
list=dir('P4\regional\GPS');
len=length(list);
disp('Step five: IONC and DCB estimate day by day!');
if license('test', 'Distrib_Computing_Toolbox')
    mypool=parpool('local',Corenum);
end
%--chose the order of spheric harmonic function
for i=3:len
    doy=list(i).name;
    total_station=length(dir(['OBS\regional\',doy,'\*.mat']));
    index=find_sp3(list_sp3,doy);
    load(['SP3\' list_sp3(index).name],'-mat');
    list_P4.gps=dir([p_opath '\GPS\',doy,'\*.mat']);
    list_P4.bds=dir([p_opath '\BDS\',doy,'\*.mat']);
    list_P4.gal=dir([p_opath '\GAL\',doy,'\*.mat']);
    list_P4.galx=dir([p_opath '\GALX\',doy,'\*.mat']);
    list_P4.glo=dir([p_opath '\GLO\',doy,'\*.mat']);
    
    if strcmp(system,'G')
            [G_R, G_S, IONC, m0, NN]= Get_nonSH_G(fig,doy,Sites_Info,sate,SDCB_REF,K,M,PG);
        fname=['M_Result\nonSH_G' doy '.mat'];
        save(fname,'G_R','G_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GC')
        [G_R, G_S, C_R, C_S, IONC,m0,NN]= Get_nonSH_GC(fig,doy,Sites_Info,sate,SDCB_REF,K,M,PG,PC);
        fname=['M_Result\nonSH_GC' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GE')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, IONC,m0,NN]= Get_nonSH_GE(fig,doy,Sites_Info,sate,SDCB_REF,K,M,PG,PE);
        fname=['M_Result\nonSH_GE' doy '.mat'];
        save(fname,'G_R','G_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GR')
        [G_R, G_S, R_R, R_S, IONC,m0,NN]= Get_nonSH_GR(fig,doy,Sites_Info,sate,SDCB_REF,K,M,PG,PR);
        fname=['M_Result\nonSH_GR' doy '.mat'];
        save(fname,'G_R','G_S','R_R','R_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GCE')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, C_R, C_S, IONC, m0, NN] = Get_nonSH_GCE(fig,doy ,Sites_Info,sate,SDCB_REF,K,M,PG,PC,PE);
        fname=['M_Result\nonSH_GCE' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GCR')
        [G_R, G_S, C_R, C_S, R_R, R_S, IONC, m0, NN] = Get_nonSH_GCR(fig,doy ,Sites_Info,sate,SDCB_REF,K,M,PG,PC,PR);
        fname=['M_Result\nonSH_GCR' doy '.mat'];
        save(fname,'G_R','G_S','C_R','C_S','R_R','R_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GER')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, R_R, R_S, IONC, m0, NN] = Get_nonSH_GER(fig,doy ,Sites_Info,sate,SDCB_REF,K,M,PG,PE,PR);
        fname=['M_Result\nonSH_GER' doy '.mat'];
        save(fname,'G_R','G_S','R_R','R_S','E_R','E_S','EX_R','EX_S','IONC','m0','NN','SDCB_REF','Sites_Info','total_station','list_P4','-mat');
    end
    if strcmp(system,'GCER')
        [G_R, G_S, E_R, E_S, EX_R, EX_S, C_R, C_S, R_R, R_S, IONC, m0, NN] = Get_nonSH_GCER(fig,doy ,Sites_Info,sate,SDCB_REF,K,M,PG,PC,PE,PR);
        fname=['M_Result\nonSH_GCER' doy '.mat'];
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