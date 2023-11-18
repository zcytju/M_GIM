function [] = read_sp3( s_ipath,s_opath,sate_mark)
%% read sp3 files to struct
%% produced from 'read_sp3.m' in M_DCB
% INPUT:
%     s_ipath: storage path of *.sp3 files
%     s_opath: storage path of satellite coordinate files
%     sate_mark: satellite status identification
%% written by Jin R et al., 2012/5/20, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/12
%% --------------------------------------------------------------------
list_obs=dir([s_ipath '/*.sp3']);
len=length(list_obs);
if len<3
 error('need at least three SP3 files !!!');
end
for i=1:len-2
    pr_obs=[s_ipath '/' list_obs(i).name];
    cu_obs=[s_ipath '/' list_obs(i+1).name];
    nx_obs=[s_ipath '/' list_obs(i+2).name];
    [gpsx1,gpsy1,gpsz1,glox1,gloy1,gloz1,galx1,galy1,galz1,bdsx1,bdsy1,bdsz1]=r_sp3(pr_obs);
    [gpsx2,gpsy2,gpsz2,glox2,gloy2,gloz2,galx2,galy2,galz2,bdsx2,bdsy2,bdsz2]=r_sp3(cu_obs);
    [gpsx3,gpsy3,gpsz3,glox3,gloy3,gloz3,galx3,galy3,galz3,bdsx3,bdsy3,bdsz3]=r_sp3(nx_obs);
    numgps=max([size(gpsx1,2),size(gpsx2,2),size(gpsx3,2)]);
    numglo=max([size(glox1,2),size(glox2,2),size(glox3,2)]);
    numgal=max([size(galx1,2),size(galx2,2),size(galx3,2)]);
    numbds=max([size(bdsx1,2),size(bdsx2,2),size(bdsx3,2)]);
    [sate.gpsx,sate.gpsy,sate.gpsz]=interplotation(numgps,gpsx1,gpsy1,gpsz1,gpsx2,gpsy2,gpsz2,gpsx3,gpsy3,gpsz3);
    [sate.glox,sate.gloy,sate.gloz]=interplotation(numglo,glox1,gloy1,gloz1,glox2,gloy2,gloz2,glox3,gloy3,gloz3);
    [sate.galx,sate.galy,sate.galz]=interplotation(numgal,galx1,galy1,galz1,galx2,galy2,galz2,galx3,galy3,galz3);
    [sate.bdsx,sate.bdsy,sate.bdsz]=interplotation(numbds,bdsx1,bdsy1,bdsz1,bdsx2,bdsy2,bdsz2,bdsx3,bdsy3,bdsz3);
    
    gpsdelete=find(sate_mark.gps==0);
    if ~isnan(gpsdelete)
        kg=length(gpsdelete);
        for j=kg:-1:1
            sate.gpsx(:,gpsdelete(j))=[];
            sate.gpsy(:,gpsdelete(j))=[];
            sate.gpsz(:,gpsdelete(j))=[];
        end
    end
    
    bdsdelete=find(sate_mark.bds==0);
    if ~isnan(bdsdelete)
        kb=length(bdsdelete);
        for j=kb:-1:1
            sate.bdsx(:,bdsdelete(j))=[];
            sate.bdsy(:,bdsdelete(j))=[];
            sate.bdsz(:,bdsdelete(j))=[];
        end
    end
    
    glodelete=find(sate_mark.glo==0);
    if ~isnan(glodelete)
        kr=length(glodelete);
        for j=kr:-1:1
            sate.glox(:,glodelete(j))=[];
            sate.gloy(:,glodelete(j))=[];
            sate.gloz(:,glodelete(j))=[];
        end
    end
    
    galdelete=find(sate_mark.gal==0);
    if ~isnan(galdelete)
        ke=length(galdelete);
        for j=ke:-1:1
            sate.galx(:,galdelete(j))=[];
            sate.galy(:,galdelete(j))=[];
            sate.galz(:,galdelete(j))=[];
        end
    end
    
    GN=list_obs(i+1).name(12:18);
    filename=strcat(GN,'sp3.mat');
    if ~isfolder(s_opath)
        mkdir(s_opath);
    end
    save([s_opath,'/',filename],'sate','-mat');
end
end
%% ----------------subfunction-----------------
function [GPSX,GPSY,GPSZ,GLOX,GLOY,GLOZ,GALX,GALY,GALZ,BDSX,BDSY,BDSZ] = r_sp3( path)
fid = fopen(path,'r');
line=fgetl(fid);
ep=0;%epoch number
day=line(12:13);
while 1
    if ~ischar(line), break, end
    line=fgetl(fid);
    if strcmp(line(1),'*') && strcmp(line(12:13),day)
        h=str2double(line(15:16));
        m=str2double(line(18:19));
        ep=h*12+round(m/5)+1;
        continue;
    end
    if strcmp(line(1),'*') && ~strcmp(line(12:13),day)
        break;
    end
    if length(line)>1 && strcmp(line(1:2),'PG')
        sv=str2double(line(3:4));
        GPSX(ep,sv) = str2double(line(5:18));
        GPSY(ep,sv) = str2double(line(19:32));
        GPSZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    if length(line)>1 && strcmp(line(1:2),'PR')
        sv=str2double(line(3:4));
        GLOX(ep,sv) = str2double(line(5:18));
        GLOY(ep,sv) = str2double(line(19:32));
        GLOZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    if length(line)>1 && strcmp(line(1:2),'PE')
        sv=str2double(line(3:4));
        GALX(ep,sv) = str2double(line(5:18));
        GALY(ep,sv) = str2double(line(19:32));
        GALZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    if length(line)>1 && strcmp(line(1:2),'PC')
        sv=str2double(line(3:4));
        BDSX(ep,sv) = str2double(line(5:18));
        BDSY(ep,sv) = str2double(line(19:32));
        BDSZ(ep,sv) = str2double(line(33:46));
        continue;
    end
end
fclose(fid);
end
%% ----------------subfunction-------------------
function [interp_x2,interp_y2,interp_z2]=interplotation(satenum,x1,y1,z1,x2,y2,z2,x3,y3,z3)
interp_x2=zeros(2880,satenum);
interp_y2=zeros(2880,satenum);
interp_z2=zeros(2880,satenum);
x2=[x1(end-3:end,:);x2;x3(1:5,:)];
y2=[y1(end-3:end,:);y2;y3(1:5,:)];
z2=[z1(end-3:end,:);z2;z3(1:5,:)];
m_t=linspace(-40,2880+40,297);
for i=1:satenum
    for j=1:288     
        tt=m_t(j:j+9);
        x=x2((j:9+j),i)';y=y2((j:9+j),i)';z=z2((j:9+j),i)';
        t0=linspace(m_t(j+4),m_t(j+5)-1,10);
        interp_x2((10*j-9):10*j,i)=interp_lag(tt,x,t0)';
        interp_y2((10*j-9):10*j,i)=interp_lag(tt,y,t0)';
        interp_z2((10*j-9):10*j,i)=interp_lag(tt,z,t0)';
    end
end
end
%% ----------------subfunction----------------
function y0 = interp_lag (x, y, x0)
n=length(x);
y0=zeros(size(x0));
for k=1:n
    t=1;
    for i=1:n
        if i~=k
            t=t.*(x0-x(i))/(x(k)-x(i));
        end
    end
    y0=y0+t*y(k);
end
end