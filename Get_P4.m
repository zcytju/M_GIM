function [] = Get_P4(path_obs,path_p4,path_sp3,Sites_Info,lim,sate_mark)
%% get smoothed P4 observations, produced from 'Get_P12.m' in M_DCB
% INPUT: 
%     path_obs: storage path of rinex files
%     path_sp3: storage path of sp3 files
%     Sites_Info: name and coordinate information of the stations
%     lim: cut-off angle
%     sate_mark: satellite status identification
%% written by Jin R et al., 2013/5/10, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/12
%% ---------------------------------------------------------------------
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
current_doy=num2str(unique(doys));
list_obs=dir([path_obs,'/',current_doy,'/*.mat']);
list_sp3=dir([path_sp3 '/*.mat']);

len=length(list_obs);
for i=1:len
    load([path_obs,'/',current_doy,'/',list_obs(i).name],'-mat');
    site=list_obs(i).name(1:end-9);
    doy=list_obs(i).name(end-8:end-4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    if sx==0 && sy==0 && sz==0
        continue;
    end %get rid of sites without receivers coordinates.
    sp3index=find_sp3(list_sp3,doy);
    load([path_sp3 '/' list_sp3(sp3index).name],'-mat');
    obs=cutobs(sate,sx,sy,sz,obs,lim,sate_mark);
    fields=fieldnames(obs);
    
    if ~isnan(find(strcmp(fields, 'GPSC1W' )))
        if ~(all(all(obs.GPSC1W==0)) || all(all(obs.GPSC2W==0)))
            [GPSP4]=GPS_prepro(obs,sate_mark);
            if all(all(GPSP4==0))
                continue;
            end
            if exist([path_p4 '/GPS/' doy],'dir')==0
                mkdir([path_p4 '/GPS/' doy]);
            end
            filenameP4=[path_p4 '/GPS/' doy '/' site doy 'P4.mat'];
            save(filenameP4,'GPSP4','-mat');
        end
    end
    
    if ~isnan(find(strcmp(fields, 'GLOC1P' )))
        if ~(all(all(obs.GLOC1P==0)) || all(all(obs.GLOC2P==0)))
            [GLOP4]=GLO_prepro(obs,sate_mark);
            if all(all(GLOP4==0))
                continue;
            end
            if exist([path_p4 '/GLO/' doy],'dir')==0
                mkdir([path_p4 '/GLO/' doy]);
            end
            filenameP4=[path_p4 '/GLO/' doy '/' site doy 'P4.mat'];
            save(filenameP4,'GLOP4','-mat');
        end
    end
    
    if ~isnan(find(strcmp(fields, 'BDSC2I' )))
        if ~(all(all(obs.BDSC2I==0)) || all(all(obs.BDSC7I==0)))
            [BDSP4]=BDS_prepro(obs,sate_mark);
            if all(all(BDSP4==0))
                continue;
            end
            if exist([path_p4 '/BDS/' doy],'dir')==0
                mkdir([path_p4 '/BDS/' doy]);
            end
            filenameP4=[path_p4 '/BDS/' doy '/' site doy 'P4.mat'];
            save(filenameP4,'BDSP4','-mat');
        end
    end
    
    if ~isnan(find(strcmp(fields, 'GALC1C' )))
        if ~(all(all(obs.GALC1C==0)) || all(all(obs.GALC5Q==0)))
            [GALP4]=GAL_prepro(obs,sate_mark);
            if all(all(GALP4==0))
                continue;
            end
            if exist([path_p4 '/GAL/' doy],'dir')==0
                mkdir([path_p4 '/GAL/' doy]);
            end
            filenameP4=[path_p4 '/GAL/' doy '/' site doy 'P4.mat'];
            save(filenameP4,'GALP4','-mat');
        end
    elseif ~isnan(find(strcmp(fields, 'GALC1X' )))
        if ~(all(all(obs.GALC1X==0)) || all(all(obs.GALC5X==0)))
            [GALXP4]=GALX_prepro(obs,sate_mark);
            if all(all(GALXP4==0))
                continue;
            end
            if exist([path_p4 '/GALX/' doy],'dir')==0
                mkdir([path_p4 '/GALX/' doy]);
            end
            filenameP4=[path_p4 '/GALX/' doy '/' site doy 'P4.mat'];
            save(filenameP4,'GALXP4','-mat');
        end
    end
    
    clear GPSP4 GLOP4 GALP4 GALXP4 BDSP4;
    clear sate;
    
end
end
%% ----------------subfunction-----------------
function obs=cutobs(sate,sx,sy,sz,obs,lim,sate_mark)
%% get the observations under specified cut-off angle
% INPUT: 
%     sate: precise coordinates of the satellites
%     sx: X coordinate of the station
%     sy: Y coordinate of the station
%     sz: Z coordinate of the station
%     obs: original observation structs
%     lim: cut-off angle
%     sate_mark: satellite status identification
% OUTPUT£º
%      obs: updated observation structs
fields=fieldnames(obs);
% cut gps obs
if ~isnan(find(strcmp(fields, 'GPSC1W' )))
    gpsline=size(obs.GPSC2W,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSL1C(gpsline+1:2880,:)=0;obs.GPSL2W(gpsline+1:2880,:)=0;
        obs.GPSC1W(gpsline+1:2880,:)=0;obs.GPSC2W(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC2W,2);
    if gpsl<gpssate
        obs.GPSL1C(:,gpsl+1:gpssate)=0;obs.GPSL2W(:,gpsl+1:gpssate)=0;
        obs.GPSC1W(:,gpsl+1:gpssate)=0;obs.GPSC2W(:,gpsl+1:gpssate)=0;
    end
    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    if gpsnum<size(obs.GPSC1W,2)
        if ~isnan(gpsdelete)
            k=length(gpsdelete);
            for i=k:-1:1
                obs.GPSC1W(:,gpsdelete(i))=[];
                obs.GPSC2W(:,gpsdelete(i))=[];
                obs.GPSL1C(:,gpsdelete(i))=[];
                obs.GPSL2W(:,gpsdelete(i))=[];
            end
        end
        obs.GPSC1W=obs.GPSC1W(:,1:gpsnum);
        obs.GPSC2W=obs.GPSC2W(:,1:gpsnum);
        obs.GPSL1C=obs.GPSL1C(:,1:gpsnum);
        obs.GPSL2W=obs.GPSL2W(:,1:gpsnum);
    end
    obs.GPSC1W(isnan(obs.GPSC1W))=0;
    obs.GPSC2W(isnan(obs.GPSC2W))=0;
    obs.GPSL1C(isnan(obs.GPSL1C))=0;
    obs.GPSL2W(isnan(obs.GPSL2W))=0;
    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;

    for i=1:gpsnum
        for j=1:2880
            if size(obs.GPSL1C, 2) ~= gpsnum
                continue;
            end
            if obs.GPSL1C(j,i)==0 || obs.GPSL2W(j,i)==0 || obs.GPSC1W(j,i)==0 || obs.GPSC2W(j,i)==0
                obs.GPSL1C(j,i)=0;obs.GPSL2W(j,i)=0;obs.GPSC1W(j,i)=0;obs.GPSC2W(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,gpsx(j,i)*1000,gpsy(j,i)*1000,gpsz(j,i)*1000);
            if el<lim
                obs.GPSC1W(j,i)=0;obs.GPSC2W(j,i)=0;obs.GPSL1C(j,i)=0;obs.GPSL2W(j,i)=0;
                continue;
            end
        end
    end
end
% cut glonass obs
if ~isnan(find(strcmp(fields, 'GLOC1P' )))
    gloline=size(obs.GLOC1P,1);
    glosate=size(sate_mark.glo,2);
    if gloline<2880
        obs.GLOL1P(gloline+1:2880,:)=0;obs.GLOL2P(gloline+1:2880,:)=0;
        obs.GLOC1P(gloline+1:2880,:)=0;obs.GLOC2P(gloline+1:2880,:)=0;
    end
    glol=size(obs.GLOC1P,2);
    if glol<glosate
        obs.GLOL1P(:,glol+1:glosate)=0;obs.GLOL2P(:,glol+1:glosate)=0;
        obs.GLOC1P(:,glol+1:glosate)=0;obs.GLOC2P(:,glol+1:glosate)=0;
    end
    
    glodelete=find(sate_mark.glo==0);
    glonum=size(sate.glox,2);    
    if glonum<size(obs.GLOC1P,2)
        if ~isnan(glodelete)
            k=length(glodelete);
            for i=k:-1:1
                obs.GLOC1P(:,glodelete(i))=[];
                obs.GLOC2P(:,glodelete(i))=[];
                obs.GLOL1P(:,glodelete(i))=[];
                obs.GLOL2P(:,glodelete(i))=[];
            end
        end
        obs.GLOC1P=obs.GLOC1P(:,1:glonum);
        obs.GLOC2P=obs.GLOC2P(:,1:glonum);
        obs.GLOL1P=obs.GLOL1P(:,1:glonum);
        obs.GLOL2P=obs.GLOL2P(:,1:glonum);
    end
    obs.GLOC1P(isnan(obs.GLOC1P))=0;
    obs.GLOC2P(isnan(obs.GLOC2P))=0;
    obs.GLOL1P(isnan(obs.GLOL1P))=0;
    obs.GLOL2P(isnan(obs.GLOL2P))=0;
    
    glox=sate.glox;gloy=sate.gloy;gloz=sate.gloz;

    for i=1:glonum
        for j=1:2880
            if size(obs.GLOL1P, 2) ~= glonum
                continue;
            end
            if obs.GLOL1P(j,i)==0 || obs.GLOL2P(j,i)==0 || obs.GLOC1P(j,i)==0 || obs.GLOC2P(j,i)==0
                obs.GLOL1P(j,i)=0;obs.GLOL2P(j,i)=0;obs.GLOC1P(j,i)=0;obs.GLOC2P(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,glox(j,i)*1000,gloy(j,i)*1000,gloz(j,i)*1000);
            if el<lim
                obs.GLOC1P(j,i)=0;obs.GLOC2P(j,i)=0;obs.GLOL1P(j,i)=0;obs.GLOL2P(j,i)=0;
                continue;
            end
        end
    end
end
% cut bds obs
if ~isnan(find(strcmp(fields, 'BDSC2I' )))
    bdsline=size(obs.BDSC2I,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSL2I(bdsline+1:2880,:)=0;obs.BDSL7I(bdsline+1:2880,:)=0;
        obs.BDSC2I(bdsline+1:2880,:)=0;obs.BDSC7I(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC2I,2);
    if bdsl<bdssate
        obs.BDSL2I(:,bdsl+1:bdssate)=0;obs.BDSL7I(:,bdsl+1:bdssate)=0;
        obs.BDSC2I(:,bdsl+1:bdssate)=0;obs.BDSC7I(:,bdsl+1:bdssate)=0;
    end
    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=sum(sate_mark.bds);
    if bdsnum<size(obs.BDSC2I,2)
        if ~isnan(bdsdelete)
            k=length(bdsdelete);
            for i=k:-1:1
                obs.BDSC2I(:,bdsdelete(i))=[];
                obs.BDSC7I(:,bdsdelete(i))=[];
                obs.BDSL2I(:,bdsdelete(i))=[];
                obs.BDSL7I(:,bdsdelete(i))=[];
            end
        end
        obs.BDSC2I=obs.BDSC2I(:,1:bdsnum);
        obs.BDSC7I=obs.BDSC7I(:,1:bdsnum);
        obs.BDSL2I=obs.BDSL2I(:,1:bdsnum);
        obs.BDSL7I=obs.BDSL7I(:,1:bdsnum);
    end
    obs.BDSC2I(isnan(obs.BDSC2I))=0;
    obs.BDSC7I(isnan(obs.BDSC7I))=0;
    obs.BDSL2I(isnan(obs.BDSL2I))=0;
    obs.BDSL7I(isnan(obs.BDSL7I))=0;
    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;

    for i=1:bdsnum
        for j=1:2880
            if size(obs.BDSL2I, 2) ~= bdsnum
                continue;
            end
            if obs.BDSL2I(j,i)==0 || obs.BDSL7I(j,i)==0 || obs.BDSC2I(j,i)==0 || obs.BDSC7I(j,i)==0
                obs.BDSL2I(j,i)=0;obs.BDSL7I(j,i)=0;obs.BDSC2I(j,i)=0;obs.BDSC7I(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,bdsx(j,i)*1000,bdsy(j,i)*1000,bdsz(j,i)*1000);
            if el<lim
                obs.BDSC2I(j,i)=0;obs.BDSC7I(j,i)=0;obs.BDSL2I(j,i)=0;obs.BDSL7I(j,i)=0;
                continue;
            end
        end
    end
end
% cut galileo obs
if ~isnan(find(strcmp(fields, 'GALC1X' )))
    galxline=size(obs.GALC1X,1);
    galsate=size(sate_mark.gal,2);
    
    if galxline<2880
        obs.GALL1X(galxline+1:2880,:)=0;obs.GALL5X(galxline+1:2880,:)=0;
        obs.GALC1X(galxline+1:2880,:)=0;obs.GALC5X(galxline+1:2880,:)=0;
    end
    gall=size(obs.GALC1X,2);
    if gall<galsate
        obs.GALL1X(:,gall+1:galsate)=0;obs.GALL5X(:,gall+1:galsate)=0;
        obs.GALC1X(:,gall+1:galsate)=0;obs.GALC5X(:,gall+1:galsate)=0;
    end
    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);
    if galnum<size(obs.GALC1X,2)
        if ~isnan(galdelete)
            k=length(galdelete);
            for i=k:-1:1
                obs.GALC1X(:,galdelete(i))=[];
                obs.GALC5X(:,galdelete(i))=[];
                obs.GALL1X(:,galdelete(i))=[];
                obs.GALL5X(:,galdelete(i))=[];
            end
        end
        obs.GALC1X=obs.GALC1X(:,1:galnum);
        obs.GALC5X=obs.GALC5X(:,1:galnum);
        obs.GALL1X=obs.GALL1X(:,1:galnum);
        obs.GALL5X=obs.GALL5X(:,1:galnum);
    end
    obs.GALC1X(isnan(obs.GALC1X))=0;
    obs.GALC5X(isnan(obs.GALC5X))=0;
    obs.GALL1X(isnan(obs.GALL1X))=0;
    obs.GALL5X(isnan(obs.GALL5X))=0;
    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;

    for i=1:galnum
        for j=1:2880
            if size(obs.GALL1X, 2) ~= galnum
                continue;
            end
            if obs.GALL1X(j,i)==0 || obs.GALL5X(j,i)==0 || obs.GALC1X(j,i)==0 || obs.GALC5X(j,i)==0
                obs.GALL1X(j,i)=0;obs.GALL5X(j,i)=0;obs.GALC1X(j,i)=0;obs.GALC5X(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,galx(j,i)*1000,galy(j,i)*1000,galz(j,i)*1000);
            if el<lim
                obs.GALC1X(j,i)=0;obs.GALC5X(j,i)=0;obs.GALL1X(j,i)=0;obs.GALL5X(j,i)=0;
                continue;
            end
        end
    end
elseif ~isnan(find(strcmp(fields, 'GALC1C' )))
    galline=size(obs.GALC1C,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALL1C(galline+1:2880,:)=0;obs.GALL5Q(galline+1:2880,:)=0;
        obs.GALC1C(galline+1:2880,:)=0;obs.GALC5Q(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC1C,2);
    if gall<galsate
        obs.GALL1C(:,gall+1:galsate)=0;obs.GALL5Q(:,gall+1:galsate)=0;
        obs.GALC1C(:,gall+1:galsate)=0;obs.GALC5Q(:,gall+1:galsate)=0;
    end
    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);
        
    if galnum<size(obs.GALC1C,2)
        if ~isnan(galdelete)
            k=length(galdelete);
            for i=k:-1:1
                obs.GALC1C(:,galdelete(i))=[];
                obs.GALC5Q(:,galdelete(i))=[];
                obs.GALL1C(:,galdelete(i))=[];
                obs.GALL5Q(:,galdelete(i))=[];
            end
        end
        obs.GALC1C=obs.GALC1C(:,1:galnum);
        obs.GALC5Q=obs.GALC5Q(:,1:galnum);
        obs.GALL1C=obs.GALL1C(:,1:galnum);
        obs.GALL5Q=obs.GALL5Q(:,1:galnum);
    end
    obs.GALC1C(isnan(obs.GALC1C))=0;
    obs.GALC5Q(isnan(obs.GALC5Q))=0;
    obs.GALL1C(isnan(obs.GALL1C))=0;
    obs.GALL5Q(isnan(obs.GALL5Q))=0;
    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;

    for i=1:galnum
        for j=1:2880
            if size(obs.GALL1C, 2) ~= galnum
                continue;
            end
            if obs.GALL1C(j,i)==0 || obs.GALL5Q(j,i)==0 || obs.GALC1C(j,i)==0 || obs.GALC5Q(j,i)==0
                obs.GALL1C(j,i)=0;obs.GALL5Q(j,i)=0;obs.GALC1C(j,i)=0;obs.GALC5Q(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,galx(j,i)*1000,galy(j,i)*1000,galz(j,i)*1000);
            if el<lim
                obs.GALC1C(j,i)=0;obs.GALC5Q(j,i)=0;obs.GALL1C(j,i)=0;obs.GALL5Q(j,i)=0;
                continue;
            end
        end
    end
end
end
%% ---------------------------subfunction----------------------------------
function [GPSPP4]=GPS_prepro(obs,sate_mark)
%% get smoothed GPS P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GPSPP4: smoothed GPS P4 observation
%% --------------------------------------------------------------------------
%____delete incomplete epoch_____
validsats=sum(sate_mark.gps);
size1=size(obs.GPSC1W,1);
size2=size(obs.GPSC1W,2);
if size2<validsats
    obs.GPSC1W(:,size2+1:validsats)=0;obs.GPSC2W(:,size2+1:validsats)=0;
    obs.GPSL1C(:,size2+1:validsats)=0;obs.GPSL2W(:,size2+1:validsats)=0;
end
GPSP4=zeros(size1,validsats);
GPSL4=zeros(size1,validsats);
GPSPP4=zeros(size1,validsats);
c=299792458;                     %__________________________speed of light
GPS_f1=1575.42*10^6;                 %_________________________________unit:Hz
GPS_f2=1227.6*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(GPS_f1-GPS_f2);       %____________________wide lane wavelength
L6=lamda_w*(obs.GPSL1C-obs.GPSL2W)-(GPS_f1*obs.GPSC1W+GPS_f2*obs.GPSC2W)/(GPS_f1+GPS_f2); %__MW observable
Li=obs.GPSL1C-GPS_f1*obs.GPSL2W/GPS_f2;
Nw=-L6;                   %_____________________wide lane ambiguity
for i=1:validsats  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GPSC1W(k,i)=0;obs.GPSC2W(k,i)=0;obs.GPSL1C(k,i)=0;obs.GPSL2W(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GPSL1C(e,i)=0;obs.GPSL2W(e,i)=0;obs.GPSC1W(e,i)=0;obs.GPSC2W(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GPSC1W(k,i)=0;obs.GPSC2W(k,i)=0;obs.GPSL1C(k,i)=0;obs.GPSL2W(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GPSC1W(k+1,i)=0;obs.GPSC2W(k+1,i)=0;
                        obs.GPSL1C(k+1,i)=0;obs.GPSL2W(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GPSC1W(l,i)=0;obs.GPSC2W(l,i)=0;
                            obs.GPSL1C(l,i)=0;obs.GPSL2W(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GPSC1W(l,i)=0;obs.GPSC2W(l,i)=0;
                            obs.GPSL1C(l,i)=0;obs.GPSL2W(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GPSC1W(k+1,i)=0;obs.GPSC2W(k+1,i)=0;
                        obs.GPSL1C(k+1,i)=0;obs.GPSL2W(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GPSC1W(l,i)=0;obs.GPSC2W(l,i)=0;
                            obs.GPSL1C(l,i)=0;obs.GPSL2W(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
    GPSP4(:,i)=obs.GPSC1W(:,i)-obs.GPSC2W(:,i);
    GPSL4(:,i)=(c/GPS_f1)*obs.GPSL1C(:,i)-(c/GPS_f2)*obs.GPSL2W(:,i);
    
    %--------------------------smoothing-------------------------
    GPSPP4(:,i)=GPSP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            GPSPP4(k,i)=GPSPP4(k,i)/t+(GPSPP4(k-1,i)+GPSL4(k-1,i)-GPSL4(k,i))*(t-1)/t;
            t=t+1;
        end
        GPSPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(GPSPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(GPSPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                GPSPP4(kk,i)=0;
            end
        end
    end
end
end

%% ------------------------subfunction----------------------------
function [GLOPP4]=GLO_prepro(obs,sate_mark)
%% get smoothed GLONASS P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GLOPP4: smoothed GLONASS P4 observation
%% ------------------------------------------------------------------
%____delete incomplete epoch_____
validsats=sum(sate_mark.glo);
size1=size(obs.GLOC1P,1);
size2=size(obs.GLOC1P,2);
if size2<validsats
    obs.GLOC1P(:,size2+1:validsats)=0;obs.GLOC2P(:,size2+1:validsats)=0;
    obs.GLOL1P(:,size2+1:validsats)=0;obs.GLOL2P(:,size2+1:validsats)=0;
end
GLOP4=zeros(size1,validsats);
GLOL4=zeros(size1,validsats);
GLOPP4=zeros(size1,validsats);
c=299792458;                     %__________________________speed of light
for i=1:validsats
    Fre=[1,-4, 5, 1, 5, 6, -2, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
    %      Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
    GLO_f1(i)=(1602+Fre(i)*0.5625)*10^6;                 %_________________________________unit:Hz
    GLO_f2(i)=(1246+Fre(i)*0.4375)*10^6;                  %_________________________________unit:Hz
    lamda_w(i)=299792458/(GLO_f1(i)-GLO_f2(i));       %____________________wide lane wavelength
    L6(:,i)=lamda_w(i)*(obs.GLOL1P(:,i)-obs.GLOL2P(:,i))-(GLO_f1(i)*obs.GLOC1P(:,i)+GLO_f2(i)*obs.GLOC2P(:,i))/(GLO_f1(i)+GLO_f2(i)); %__MW observable
    Li(:,i)=obs.GLOL1P(:,i)-GLO_f1(i)*obs.GLOL2P(:,i)/GLO_f2(i);
    Nw(:,i)=-L6(:,i);               %_____________________wide lane ambiguity
end
for i=1:validsats  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GLOC1P(k,i)=0;obs.GLOC2P(k,i)=0;obs.GLOL1P(k,i)=0;obs.GLOL2P(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GLOL1P(e,i)=0;obs.GLOL2P(e,i)=0;obs.GLOC1P(e,i)=0;obs.GLOC2P(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GLOC1P(k,i)=0;obs.GLOC2P(k,i)=0;obs.GLOL1P(k,i)=0;obs.GLOL2P(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GLOC1P(k+1,i)=0;obs.GLOC2P(k+1,i)=0;
                        obs.GLOL1P(k+1,i)=0;obs.GLOL2P(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GLOC1P(l,i)=0;obs.GLOC2P(l,i)=0;
                            obs.GLOL1P(l,i)=0;obs.GLOL2P(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GLOC1P(l,i)=0;obs.GLOC2P(l,i)=0;
                            obs.GLOL1P(l,i)=0;obs.GLOL2P(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GLOC1P(k+1,i)=0;obs.GLOC2P(k+1,i)=0;
                        obs.GLOL1P(k+1,i)=0;obs.GLOL2P(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GLOC1P(l,i)=0;obs.GLOC2P(l,i)=0;
                            obs.GLOL1P(l,i)=0;obs.GLOL2P(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
    GLOP4(:,i)=obs.GLOC1P(:,i)-obs.GLOC2P(:,i);
    GLOL4(:,i)=(c/GLO_f1(i))*obs.GLOL1P(:,i)-(c/GLO_f2(i))*obs.GLOL2P(:,i);
    
    %--------------------------smoothing-------------------------
    GLOPP4(:,i)=GLOP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            GLOPP4(k,i)=GLOPP4(k,i)/t+(GLOPP4(k-1,i)+GLOL4(k-1,i)-GLOL4(k,i))*(t-1)/t;
            t=t+1;
        end
        GLOPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(GLOPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(GLOPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                GLOPP4(kk,i)=0;
            end
        end
    end
end
end

%% ---------------------------subfunction-------------------------
function [BDSPP4]=BDS_prepro(obs,sate_mark)
%% get smoothed BDS P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     BDSPP4: smoothed BDS P4 observation
%% -----------------------------------------------------------------
%____delete incomplete epoch_____
validsats=sum(sate_mark.bds);
size1=size(obs.BDSC2I,1);
size2=size(obs.BDSC2I,2);
if size2<validsats
    obs.BDSC2I(:,size2+1:validsats)=0;obs.BDSC7I(:,size2+1:validsats)=0;
    obs.BDSL2I(:,size2+1:validsats)=0;obs.BDSL7I(:,size2+1:validsats)=0;
end
BDSP4=zeros(size1,validsats);
BDSL4=zeros(size1,validsats);
BDSPP4=zeros(size1,validsats);
c=299792458;                     %__________________________speed of light
BDS_f2=1561.098*10^6;                 %_________________________________unit:Hz
BDS_f7=1207.140*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(BDS_f2-BDS_f7);       %____________________wide lane wavelength
L6=lamda_w*(obs.BDSL2I-obs.BDSL7I)-(BDS_f2*obs.BDSC2I+BDS_f7*obs.BDSC7I)/(BDS_f2+BDS_f7); %__MW observable
Li=obs.BDSL2I-BDS_f2*obs.BDSL7I/BDS_f7;
Nw=-L6;                   %_____________________wide lane ambiguity
for i=1:validsats  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.BDSC2I(k,i)=0;obs.BDSC7I(k,i)=0;obs.BDSL2I(k,i)=0;obs.BDSL7I(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.BDSL2I(e,i)=0;obs.BDSL7I(e,i)=0;obs.BDSC2I(e,i)=0;obs.BDSC7I(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.BDSC2I(k,i)=0;obs.BDSC7I(k,i)=0;obs.BDSL2I(k,i)=0;obs.BDSL7I(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.BDSC2I(k+1,i)=0;obs.BDSC7I(k+1,i)=0;
                        obs.BDSL2I(k+1,i)=0;obs.BDSL7I(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.BDSC2I(l,i)=0;obs.BDSC7I(l,i)=0;
                            obs.BDSL2I(l,i)=0;obs.BDSL7I(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.BDSC2I(l,i)=0;obs.BDSC7I(l,i)=0;
                            obs.BDSL2I(l,i)=0;obs.BDSL7I(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.BDSC2I(k+1,i)=0;obs.BDSC7I(k+1,i)=0;
                        obs.BDSL2I(k+1,i)=0;obs.BDSL7I(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.BDSC2I(l,i)=0;obs.BDSC7I(l,i)=0;
                            obs.BDSL2I(l,i)=0;obs.BDSL7I(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
    BDSP4(:,i)=obs.BDSC2I(:,i)-obs.BDSC7I(:,i);
    BDSL4(:,i)=(c/BDS_f2)*obs.BDSL2I(:,i)-(c/BDS_f7)*obs.BDSL7I(:,i);
    
    %--------------------------smoothing-------------------------
    BDSPP4(:,i)=BDSP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            BDSPP4(k,i)=BDSPP4(k,i)/t+(BDSPP4(k-1,i)+BDSL4(k-1,i)-BDSL4(k,i))*(t-1)/t;
            t=t+1;
        end
        BDSPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(BDSPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(BDSPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                BDSPP4(kk,i)=0;
            end
        end
    end
end
end

%% --------------------------subfunction--------------------------
function [GALPP4]=GAL_prepro(obs,sate_mark)
%% get smoothed Galileo P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GALPP4: smoothed Galileo P4 observation
%% ------------------------------------------------------------------
%____delete incomplete epoch_____
validsats=sum(sate_mark.gal);
size1=size(obs.GALC1C,1);
size2=size(obs.GALC1C,2);
if size2<validsats
    obs.GALL1C(:,size2+1:validsats)=0;obs.GALL5Q(:,size2+1:validsats)=0;
    obs.GALC1C(:,size2+1:validsats)=0;obs.GALC5Q(:,size2+1:validsats)=0;
end
GALP4=zeros(size1,validsats);
GALL4=zeros(size1,validsats);
GALPP4=zeros(size1,validsats);
c=299792458;                     %__________________________speed of light
GAL_f1=1575.42*10^6;                 %_________________________________unit:Hz
GAL_f5=1176.45*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(GAL_f1-GAL_f5);       %____________________wide lane wavelength
L6=lamda_w*(obs.GALL1C-obs.GALL5Q)-(GAL_f1*obs.GALC1C+GAL_f5*obs.GALC5Q)/(GAL_f1+GAL_f5); %__MW observable
Li=obs.GALL1C-GAL_f1*obs.GALL5Q/GAL_f5;
Nw=-L6;                   %_____________________wide lane ambiguity
for i=1:validsats  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1C(k,i)=0;obs.GALC5Q(k,i)=0;obs.GALL1C(k,i)=0;obs.GALL5Q(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GALL1C(e,i)=0;obs.GALL5Q(e,i)=0;obs.GALC1C(e,i)=0;obs.GALC5Q(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1C(k,i)=0;obs.GALC5Q(k,i)=0;obs.GALL1C(k,i)=0;obs.GALL5Q(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1C(k+1,i)=0;obs.GALC5Q(k+1,i)=0;
                        obs.GALL1C(k+1,i)=0;obs.GALL5Q(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GALC1C(l,i)=0;obs.GALC5Q(l,i)=0;
                            obs.GALL1C(l,i)=0;obs.GALL5Q(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GALC1C(l,i)=0;obs.GALC5Q(l,i)=0;
                            obs.GALL1C(l,i)=0;obs.GALL5Q(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1C(k+1,i)=0;obs.GALC5Q(k+1,i)=0;
                        obs.GALL1C(k+1,i)=0;obs.GALL5Q(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GALC1C(l,i)=0;obs.GALC5Q(l,i)=0;
                            obs.GALL1C(l,i)=0;obs.GALL5Q(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
    
    GALP4(:,i)=obs.GALC1C(:,i)-obs.GALC5Q(:,i);
    GALL4(:,i)=(c/GAL_f1)*obs.GALL1C(:,i)-(c/GAL_f5)*obs.GALL5Q(:,i);
    
    %--------------------------smoothing-------------------------
    GALPP4(:,i)=GALP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            GALPP4(k,i)=GALPP4(k,i)/t+(GALPP4(k-1,i)+GALL4(k-1,i)-GALL4(k,i))*(t-1)/t;
            t=t+1;
        end
        GALPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(GALPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(GALPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                GALPP4(kk,i)=0;
            end
        end
    end
end
end

%% --------------------------subfunction----------------------------
function [GALXPP4]=GALX_prepro(obs,sate_mark)
%% get smoothed Galileo P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GALXPP4: smoothed Galileo P4 observation
%% --------------------------------------------------------------------
%____delete incomplete epoch_____
validsats=sum(sate_mark.gal);
size1=size(obs.GALC1X,1);
size2=size(obs.GALC1X,2);
if size2<validsats
    obs.GALL1X(:,size2+1:validsats)=0;obs.GALL5X(:,size2+1:validsats)=0;
    obs.GALC1X(:,size2+1:validsats)=0;obs.GALC5X(:,size2+1:validsats)=0;
end
GALXP4=zeros(size1,validsats);
GALXL4=zeros(size1,validsats);
GALXPP4=zeros(size1,validsats);
c=299792458;                     %__________________________speed of light
GALX_f1=1575.42*10^6;                 %_________________________________unit:Hz
GALX_f5=1176.45*10^6;                 %_________________________________unit:Hz
lamda_w=299792458/(GALX_f1-GALX_f5);       %____________________wide lane wavelength
L6=lamda_w*(obs.GALL1X-obs.GALL5X)-(GALX_f1*obs.GALC1X+GALX_f5*obs.GALC5X)/(GALX_f1+GALX_f5); %__MW observable
Li=obs.GALL1X-GALX_f1*obs.GALL5X/GALX_f5;
Nw=-L6;                   %_____________________wide lane ambiguity
for i=1:validsats  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1X(k,i)=0;obs.GALC5X(k,i)=0;obs.GALL1X(k,i)=0;obs.GALL5X(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GALL1X(e,i)=0;obs.GALL5X(e,i)=0;obs.GALC1X(e,i)=0;obs.GALC5X(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1X(k,i)=0;obs.GALC5X(k,i)=0;obs.GALL1X(k,i)=0;obs.GALL5X(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1X(k+1,i)=0;obs.GALC5X(k+1,i)=0;
                        obs.GALL1X(k+1,i)=0;obs.GALL5X(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GALC1X(l,i)=0;obs.GALC5X(l,i)=0;
                            obs.GALL1X(l,i)=0;obs.GALL5X(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GALC1X(l,i)=0;obs.GALC5X(l,i)=0;
                            obs.GALL1X(l,i)=0;obs.GALL5X(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1X(k+1,i)=0;obs.GALC5X(k+1,i)=0;
                        obs.GALL1X(k+1,i)=0;obs.GALL5X(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GALC1X(l,i)=0;obs.GALC5X(l,i)=0;
                            obs.GALL1X(l,i)=0;obs.GALL5X(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
    GALXP4(:,i)=obs.GALC1X(:,i)-obs.GALC5X(:,i);
    GALXL4(:,i)=(c/GALX_f1)*obs.GALL1X(:,i)-(c/GALX_f5)*obs.GALL5X(:,i);
    
    %--------------------------smoothing-------------------------
    GALXPP4(:,i)=GALXP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            GALXPP4(k,i)=GALXPP4(k,i)/t+(GALXPP4(k-1,i)+GALXL4(k-1,i)-GALXL4(k,i))*(t-1)/t;
            t=t+1;
        end
        GALXPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(GALXPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(GALXPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                GALXPP4(kk,i)=0;
            end
        end
    end
end
end
%% 
function arc = Get_arc(array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
len=length(array);
arc=[];
for i=1:len
         if i==len
             if array(i)~=0
                 arc=[arc,i];
             end
             continue;
         end
         if i==1&&array(i)~=0
             arc=[arc,i];
         end
         if array(i)==0&&array(i+1)~=0
             arc=[arc,i+1];
             continue;
         end
         if array(i)~=0&&array(i+1)==0
             arc=[arc,i];
             continue;
         end
end
if len==0,return,end
arc=reshape(arc,2,[]);
arc=arc';
end