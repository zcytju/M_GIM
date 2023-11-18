%----------------------DCB estimate based on multi stations----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

%-----------------------------Setting-------------------------------------- 
%--rinex
r_ipath=input('Please input Rinex files path:','s');
%--sp3
s_ipath=input('Please input SP3 files path:','s');
%--ionex
i_ipath=input('Please input ionex files path:','s');
%--el
lim=str2double(input('Please input elevation angle threshold(unit:degree)(10 degree is recommended):','s'));
%--order
order=str2double(input('Please input the order of spheric harmonic function (4 order is recommended):','s'));
%Step one--------------------Read RINEX files------------------------------
disp('MDCB(multi-stations) starts running!')
%--L1,L2,P1,P2 and coordinates of GPS receviers are obtained by this step
r_opath='M_OBS';    
disp('Step one: read rinex files !')
Sites_Info=read_rinex(r_ipath,r_opath);
disp('Step one: completing !')
%Step two--------------------Read SP3 files--------------------------------
%--Get satellites coordinates and interpolation
%--SP3 files of the previous day and the next day are needed for interpolation
s_opath='SP3';
disp('Step two: read SP3 files !')
read_sp3(s_ipath,s_opath);
disp('Step two: completing !')
%Step three------------------Read ionex files------------------------------
%--Get reference information
disp('Step three: read ionex files !')
[SDCB_REF,Sites_Info]= read_ionex(i_ipath,Sites_Info);
save('SDCB_REF.mat','SDCB_REF');
save('Sites_Info.mat','Sites_Info');
disp('Step three: completing !')
% 
%Step four------------------Ionosphere Observations-----------------------
%--get smoothed P4 observations
r_opath='M_OBS';
s_opath='SP3';
disp('Step four: data pre-processing and get ionosphere observations !')
Get_P12(r_opath,s_opath,Sites_Info,lim*pi/180,0);
disp('Step four: completing !')
%Step five-------------------Get M_DCB-------------------------------------
list_sp3=dir('SP3/*.mat');
list=dir('M_P4');
len=length(list);
disp('Step five: DCB estimate day by day!')
%--chose the order of spheric harmonic function
for i=3:len
    doy=list(i).name;
    index=find_sp3(list_sp3,doy);
    load(['SP3/' list_sp3(index).name],'-mat');
    [DCB_R DCB_S IONC]= Get_MDCB(doy,Sites_Info,sate,SDCB_REF,order);
    if exist('M_Result','dir')==0
        mkdir('M_Result');
    end
    fname=['M_Result/M_R' doy '.mat'];
    save(fname,'DCB_R','DCB_S','IONC','-mat');
    clear sate; 
    disp(['DCB estimate of doy ' doy ' complete!']);
end
disp('Step five: completing!')
