%----------------------DCB estimate based on single stations----------------
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
disp('MDCB(single-stations) starts running!')
%--L1,L2,P1,P2 and coordinates of GPS receviers are obtained by this step
r_opath='S_OBS'; 
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
%Step four------------------Ionosphere Observations-----------------------
%--get smoothed P4 observations
disp('Step four: data pre-processing and get ionosphere observations !')
Get_P12(r_opath,s_opath,Sites_Info,lim*pi/180,1);
disp('Step four: completing !')
% %Step five-------------------Get M_DCB-------------------------------------
list_sp3=dir('SP3/*.mat');
list=dir('S_P4');
len=length(list);
if exist('S_Result','dir')==0
    mkdir('S_Result');
end
disp('Step five: DCB estimate site by site and day by day!')
for i=3:len
    site=list(i).name;
    if exist(['S_Result/' site],'dir')==0
        mkdir(['S_Result/' site]);
    end
    list_obs=dir(['S_P4/' site '/*.mat']);
    n_d=length(list_obs);
    for j=1:n_d
        doy=list_obs(j).name(5:9);
        [DCB_R DCB_S IONC]= Get_SDCB(site,doy,j,Sites_Info,s_opath,SDCB_REF,order);
        fname=['S_Result/' site '/S' doy '.mat'];
        save(fname,'DCB_R','DCB_S','IONC','-mat');
        disp([ 'DCB estimate of GPS station ' site '  doy ' doy ' complete!']);
    end
end
disp('Step five: completing!')
