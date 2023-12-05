function [G_R, G_S, E_R, E_S, EX_R, EX_S, IONC, m0, NN] = Get_SH_GE(fig,doy ,Sites_Info,sate,SDCB_REF,order,PG,PE,sate_mark)
%%  estimate satellite and receiver DCBs, ionospheric parameters, et al.
%%  produced from 'Get_MDCB.m' in M_DCB 
% INPUT:
%     fig: group number of products
%     doy: year and doy of year
%     Sites_Info: name and coordinate information of the stations
%     sate: precise coordinates of the satellites
%     SDCB_REF: reference satellite DCBs
%     order: the order of SH model
%     PG,PE: weight of different system observations
% OUTPUT:
%     G_R, G_S, E_R, E_S, EX_R, EX_S: estimated receiver and satellite DCBs
%     IONC: ionospheric parameters
%     m0: standard deviation
%     NN: covariance matrices
%% written by Jin R et al., 2012/5/26, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/14
%% --------------------------------------------------------------------------
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
%% check the gps data
gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
path_G=['P4/global/GPS/' doy];
list_gps=dir([path_G '/*.mat']);
G_n_r=length(list_gps);%the number of receivers
%--check the number of each satellite's observations 
gpsnum=sum(sate_mark.gps);
G_PRN=linspace(0,0,gpsnum);
G_S=linspace(0,0,gpsnum);
for i=1:G_n_r
    load([path_G '/' list_gps(i).name],'-mat');
    for j=1:gpsnum
        for k=1:2880
            if GPSP4(k,j)~=0
                G_PRN(j)=G_PRN(j)+1;
            end
        end
    end
    clear GPSP4;
end
gps_d_sat=find(G_PRN==0);
if isempty(gps_d_sat)
    G_n_s=gpsnum;
else
    G_n_s=gpsnum-length(gps_d_sat);%the number of satellites
    disp(['doy ', doy ,' GPS PRN ',num2str(gps_d_sat) ,' have no observations.']);
    for k=length(gps_d_sat):-1:1
        gpsx(:,gps_d_sat(k))=[];gpsy(:,gps_d_sat(k))=[];gpsz(:,gps_d_sat(k))=[];
    end
end

if G_n_s==gpsnum
    G_Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy); 
    G_Wx=-sum(SDCB_REF.gps(index,gps_d_sat));
    G_S(gps_d_sat)=SDCB_REF.gps(index,gps_d_sat);
end

%% check the galileo data
galx=sate.galx;galy=sate.galy;galz=sate.galz;
path_E=['P4/global/GAL/' doy];
list_gal=dir([path_E '/*.mat']);
E_n_r=length(list_gal);%the number of receivers
%--check the number of each satellite's observations 
galnum=sum(sate_mark.gal);
E_PRN=linspace(0,0,galnum);
E_S=linspace(0,0,galnum);
for i=1:E_n_r
    load([path_E '/' list_gal(i).name],'-mat');
    for j=1:galnum
        for k=1:2880
            if GALP4(k,j)~=0
                E_PRN(j)=E_PRN(j)+1;
            end
        end
    end
    clear GALP4;
end
gal_d_sat=find(E_PRN==0,1);
if isempty(gal_d_sat)
    E_n_s=galnum;
else
    E_n_s=galnum-length(gal_d_sat);%the number of satellites
    disp(['doy ', doy ,' GAL PRN ',num2str(gal_d_sat) ,' have no observations.']);
    for k=length(gal_d_sat):-1:1
        galx(:,gal_d_sat(k))=[];galy(:,gal_d_sat(k))=[];galz(:,gal_d_sat(k))=[];
    end
end
if E_n_s==galnum
    E_Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy); 
    E_Wx=-sum(SDCB_REF.gal(index,gal_d_sat));
    E_S(gal_d_sat)=SDCB_REF.gal(index,gal_d_sat);
end

%% check the galileo(x) data
galxx=sate.galx;galxy=sate.galy;galxz=sate.galz;
path_EX=['P4/global/GALX/' doy];
list_galx=dir([path_EX '/*.mat']);
EX_n_r=length(list_galx);%the number of receivers
%--check the number of each satellite's observations 
EX_PRN=linspace(0,0,galnum);
EX_S=linspace(0,0,galnum);
for i=1:EX_n_r
    load([path_EX '/' list_galx(i).name],'-mat');
    for j=1:galnum
        for k=1:2880
            if GALXP4(k,j)~=0
                EX_PRN(j)=EX_PRN(j)+1;
            end
        end
    end
    clear GALXP4;
end
galx_d_sat=find(EX_PRN==0,1);
if isempty(galx_d_sat)
    EX_n_s=galnum;
else
    EX_n_s=galnum-length(galx_d_sat);%the number of satellites
    disp(['doy ', doy ,' GALX PRN ',num2str(galx_d_sat) ,' have no observations.']);
    for k=length(galx_d_sat):-1:1
        galxx(:,galx_d_sat(k))=[];galxy(:,galx_d_sat(k))=[];galxz(:,galx_d_sat(k))=[];
    end
end
if EX_n_s==galnum
    EX_Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy); 
    EX_Wx=-sum(SDCB_REF.galx(index,galx_d_sat));
    EX_S(galx_d_sat)=SDCB_REF.galx(index,galx_d_sat);
end

%% --chose the order of spheric harmonic function
%order=str2double(input('Please input the order of spheric harmonic function (4 order is recommended):','s'));
%--LS estimate
num=(order+1)^2*fig+G_n_s+G_n_r+E_n_s+E_n_r+EX_n_s+EX_n_r;
N=zeros(num,num);
U=zeros(num,1);
L=0; sizel=0;
C_GPS=linspace(0,0,num);
C_GAL=linspace(0,0,num);
C_GALX=linspace(0,0,num);
C_GPS(G_n_r+1:G_n_r+G_n_s)=ones(1,G_n_s);
C_GAL(G_n_r+G_n_s+E_n_r+1:G_n_r+G_n_s+E_n_r+E_n_s)=ones(1,E_n_s);
C_GALX(G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+1:G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+EX_n_s)=ones(1,EX_n_s);

for i=1:G_n_r
    load([path_G '/' list_gps(i).name],'-mat');
    if ~isempty(gps_d_sat)
        for k=length(gps_d_sat):-1:1
            GPSP4(:,gps_d_sat(k))=[];
        end
    end
    site=list_gps(i).name(1:4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    [sN,sl]=Get_GPSMatrix(fig,GPSP4,gpsx,gpsy,gpsz,sx,sy,sz,G_n_r,E_n_r,EX_n_r,G_n_s,E_n_s,EX_n_s,i,order);
    N=N+sN'*sN*PG;
    U=U+sN'*sl*PG;
    sizel=sizel+length(sl);
    L=L+sl'*sl*PG;
    clear GPSP4;
    disp(['1.----- [ ',num2str(i),' / ',num2str(G_n_r),' ] ',num2str(i/G_n_r*100),'% GPS data has constructed !']);
end

for i=1:E_n_r
    load([path_E '/' list_gal(i).name],'-mat');
    if ~isempty(gal_d_sat)
        for k=length(gal_d_sat):-1:1
            GALP4(:,gal_d_sat(k))=[];
        end
    end
    site=list_gal(i).name(1:4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    [sN,sl]=Get_GALMatrix(fig,GALP4,galx,galy,galz,sx,sy,sz,G_n_r,E_n_r,EX_n_r,G_n_s,E_n_s,EX_n_s,i,order);
    N=N+sN'*sN*PE;
    U=U+sN'*sl*PE;
    sizel=sizel+length(sl);
    L=L+sl'*sl*PE;
    clear GALP4;
    disp(['2.----- [ ',num2str(i),' / ',num2str(E_n_r+EX_n_r),' ] ',num2str(i/(E_n_r+EX_n_r)*100),'% Galileo data has constructed !']);
end

for i=1:EX_n_r
    load([path_EX '/' list_galx(i).name],'-mat');
    if ~isempty(galx_d_sat)
        for k=length(galx_d_sat):-1:1
            GALXP4(:,galx_d_sat(k))=[];
        end
    end
    site=list_galx(i).name(1:4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    [sN,sl]=Get_GALXMatrix(fig,GALXP4,galxx,galxy,galxz,sx,sy,sz,G_n_r,E_n_r,EX_n_r,G_n_s,E_n_s,EX_n_s,i,order);
    N=N+sN'*sN*PE;
    U=U+sN'*sl*PE;
    sizel=sizel+length(sl);
    L=L+sl'*sl*PE;
    clear GALXP4;
    disp(['2.----- [ ',num2str(E_n_r+i),' / ',num2str(E_n_r+EX_n_r),' ] ',num2str((E_n_r+i)/(E_n_r+EX_n_r)*100),'% Galileo data has constructed !']);
end

N=N+C_GPS'*C_GPS+C_GAL'*C_GAL+C_GALX'*C_GALX;
U=U+C_GPS'*G_Wx+C_GAL'*E_Wx+C_GALX'*EX_Wx;
L=L+G_Wx'*G_Wx+E_Wx'*E_Wx+EX_Wx'*EX_Wx;
R=pinv(N)*U;
G_R=R(1:G_n_r)*10^9/299792458;
temp_gps=linspace(1,gpsnum,gpsnum);
temp_gps(gps_d_sat)=[];
G_S(temp_gps)=R(G_n_r+1:G_n_r+G_n_s)*10^9/299792458;

E_R=R(G_n_r+G_n_s+1:G_n_r+G_n_s+E_n_r)*10^9/299792458;
temp_gal=linspace(1,galnum,galnum);
temp_gal(gal_d_sat)=[];
E_S(temp_gal)=R(G_n_r+G_n_s+E_n_r+1:G_n_r+G_n_s+E_n_r+E_n_s)*10^9/299792458;

EX_R=R(G_n_r+G_n_s+E_n_r+E_n_s+1:G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r)*10^9/299792458;
temp_galx=linspace(1,galnum,galnum);
temp_galx(galx_d_sat)=[];
EX_S(temp_galx)=R(G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+1:G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+EX_n_s)*10^9/299792458;

IONC=R(G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+EX_n_s+1:end);

V=L-R'*U;
f=sizel-num;
m0=sqrt(V/f);
NN=N(G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+EX_n_s+1:end,G_n_r+G_n_s+E_n_r+E_n_s+EX_n_r+EX_n_s+1:end);
end

%% ------------------------------sub_function--------------------------------
function [M,l]=Get_GPSMatrix(fig,GPSP4,x,y,z,sx,sy,sz,gps_n_r,gal_n_r,galx_n_r,gps_n_s,gal_n_s,galx_n_s,ith,order)
M=[];
l=[];
num=(order+1)^2;
[sb,sl]=XYZtoBLH(sx,sy,sz);
figt=2880/fig;
for i=1:fig
    for j=1:gps_n_s                %----j is satellite number
        parfor k=figt*i-(figt-1):figt*i %----k is epoch number
            if GPSP4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,num*fig+gps_n_s+gps_n_r+gal_n_s+gal_n_r+galx_n_s+galx_n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j)*1000,y(k,j)*1000,z(k,j)*1000); %----sx,sy,sz station coordinate ; x,y,z satellite coordinate
            IPPz=asin(6371000*sin(pi/2-E)/(6371000+450000));%-----SLM 
            t_r=30*(k-1)*pi/43200;
            [b,s]=Get_IPP(E,A,sb,sl,IPPz,t_r);%
            M_col(ith)=(-9.52437)*cos(IPPz);   %----station dcb coefficient
            M_col(gps_n_r+j)=(-9.52437)*cos(IPPz); %----satallite dcb coefficient
            st=num*(i-1)+gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+galx_n_s+1;
            ed=num*i+gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+galx_n_s;
            M_col(st:ed)=Get_SH(b,s,order); %----spherical harmonic model
            M_scol=sparse(M_col);
            M=[M;M_scol]; 
            l=[l;GPSP4(k,j)*(-9.52437)*cos(IPPz)];
        end
    end
end
end
%% ------------------------------sub_function--------------------------------
function [M,l]=Get_GALMatrix(fig,GALP4,x,y,z,sx,sy,sz,gps_n_r,gal_n_r,galx_n_r,gps_n_s,gal_n_s,galx_n_s,ith,order)
M=[];
l=[];
num=(order+1)^2;
[sb,sl]=XYZtoBLH(sx,sy,sz);
figt=2880/fig;
for i=1:fig
    for j=1:gal_n_s                %-----------------------j is satellite number
        parfor k=figt*i-(figt-1):figt*i %-------------------------k is epoch number
            if GALP4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,num*fig+gps_n_s+gps_n_r+gal_n_s+gal_n_r+galx_n_s+galx_n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j)*1000,y(k,j)*1000,z(k,j)*1000); %----sx,sy,sz station coordinate ; x,y,z satellite coordinate
            IPPz=asin(6371000*sin(pi/2-E)/(6371000+450000));%-----SLM 
            t_r=30*(k-1)*pi/43200;
            [b,s]=Get_IPP(E,A,sb,sl,IPPz,t_r);%
            M_col(gps_n_r+gps_n_s+ith)=(-7.76751393)*cos(IPPz);   %-----station dcb coefficient
            M_col(gps_n_r+gps_n_s+gal_n_r+j)=(-7.76751393)*cos(IPPz); %---satallite dcb coefficient
            st=num*(i-1)+gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+galx_n_s+1;
            ed=num*i+gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+galx_n_s;
            M_col(st:ed)=Get_SH(b,s,order); %----spherical harmonic model 
            M_scol=sparse(M_col);
            M=[M;M_scol]; 
            l=[l;GALP4(k,j)*(-7.76751393)*cos(IPPz)];
        end
    end
end
end
%% ------------------------------sub_function--------------------------------
function [M,l]=Get_GALXMatrix(fig,GALXP4,x,y,z,sx,sy,sz,gps_n_r,gal_n_r,galx_n_r,gps_n_s,gal_n_s,galx_n_s,ith,order)
M=[];
l=[];
num=(order+1)^2;
[sb,sl]=XYZtoBLH(sx,sy,sz);
figt=2880/fig;
for i=1:fig
    for j=1:galx_n_s                %----j is satellite number
        parfor k=figt*i-(figt-1):figt*i %----k is epoch number
            if GALXP4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,num*fig+gps_n_s+gps_n_r+gal_n_s+gal_n_r+galx_n_s+galx_n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j)*1000,y(k,j)*1000,z(k,j)*1000); %----sx,sy,sz station coordinate ; x,y,z satellite coordinate
            IPPz=asin(6371000*sin(pi/2-E)/(6371000+450000));%-----SLM 
            t_r=30*(k-1)*pi/43200;
            [b,s]=Get_IPP(E,A,sb,sl,IPPz,t_r);%
            M_col(gps_n_r+gps_n_s+gal_n_r+gal_n_s+ith)=(-7.76751393)*cos(IPPz);   %----station dcb coefficient
            M_col(gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+j)=(-7.76751393)*cos(IPPz); %----satallite dcb coefficient
            st=num*(i-1)+gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+galx_n_s+1;
            ed=num*i+gps_n_r+gps_n_s+gal_n_r+gal_n_s+galx_n_r+galx_n_s;
            M_col(st:ed)=Get_SH(b,s,order); %----spherical harmonic model
            M_scol=sparse(M_col);
            M=[M;M_scol]; 
            l=[l;GALXP4(k,j)*(-7.76751393)*cos(IPPz)];
        end
    end
end
end
%% -----------------------------------------------------------------------------
function cof_P=Get_SH(b,s,order)
cof_P=linspace(0,0,(order+1)^2);
ms=linspace(s,order*s,order);
i=1;
x=sin(b);
for n=0:order
    P=legendre(n,x);
    for m=0:n
        if m==0
            cof_P(i)=P(m+1)*norm(n,m);                    %------------an0
        else
            cof_P(i)=P(m+1)*norm(n,m)*cos(ms(m));         %------------anm
            i=i+1;
            cof_P(i)=P(m+1)*norm(n,m)*sin(ms(m));         %------------bnm
        end
        i=i+1;
    end
end
end
%------------------------------sub_function--------------------------------
function N=norm(n,m)
if m==0
    N=sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
else
    N=sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
end
end
