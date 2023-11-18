function [] = read_sp3( s_ipath,s_opath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
list_obs=dir([s_ipath '/*.sp3']);
len=length(list_obs);
if len<3
 error('need at least three SP3 files; see Readme');
end
for i=1:len-2;
    pr_obs=[s_ipath '/' list_obs(i).name];
    cu_obs=[s_ipath '/' list_obs(i+1).name];
    nx_obs=[s_ipath '/' list_obs(i+2).name];
    [x1,y1,z1]=r_sp3(pr_obs);
    [x2,y2,z2]=r_sp3(cu_obs);
    [x3,y3,z3]=r_sp3(nx_obs);
    [sate.x,sate.y,sate.z]=interplotation(x1,y1,z1,x2,y2,z2,x3,y3,z3);
    GN=str2double(list_obs(i+1).name(4:7));
    Day=str2double(list_obs(i+1).name(8));
    [y,doy]=GwToDoy(GN,Day); 
    filename=strcat(num2str(y,'%04d'),num2str(doy,'%03d'),'sp3.mat');
    if exist(s_opath,'dir')==0
        mkdir(s_opath);
    end
    save([s_opath,'/',filename],'sate','-mat');
end
end
%----------------subfunction-----------------------------------------------
function [X,Y,Z] = r_sp3( path)
%READ_SP3 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(path,'r');
line=fgetl(fid);
ep=0;%epoch number
X=zeros(96,32);
Y=zeros(96,32);
Z=zeros(96,32);
while 1
    if ~ischar(line), break, end;
    line=fgetl(fid);
    if strcmp(line(1),'*')
        h=str2double(line(15:16));
        m=str2double(line(18:19));
        ep=h*4+round(m/15)+1;
        continue;
    end        %---------------------get epoch number
    if length(line)>1&&strcmp(line(1:2),'PG')
        sv=str2double(line(3:4));
        X(ep,sv) = str2double(line(5:18));
        Y(ep,sv) = str2double(line(19:32));
        Z(ep,sv) = str2double(line(33:46));
        continue;
    end
end
fclose(fid);
end
%----------------subfunction-----------------------------------------------
function [interp_x2,interp_y2,interp_z2]=interplotation(x1,y1,z1,x2,y2,z2,x3,y3,z3)
interp_x2=zeros(2880,32);
interp_y2=zeros(2880,32);
interp_z2=zeros(2880,32);
x2=[x1(93:96,:);x2;x3(1:5,:)];
y2=[y1(93:96,:);y2;y3(1:5,:)];
z2=[z1(93:96,:);z2;z3(1:5,:)];
m_t=linspace(-120,3000,105);
for i=1:32;
    for j=1:96     
        tt=m_t(j:j+9);
        x=x2((j:9+j),i)';y=y2((j:9+j),i)';z=z2((j:9+j),i)';
        t0=linspace(m_t(j+4),m_t(j+5)-1,30);
        interp_x2((30*j-29):30*j,i)=interp_lag(tt,x,t0)';
        interp_y2((30*j-29):30*j,i)=interp_lag(tt,y,t0)';
        interp_z2((30*j-29):30*j,i)=interp_lag(tt,z,t0)';
    end
end
end
%----------------subfunction-----------------------------------------------
function y0 = interp_lag (x, y, x0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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