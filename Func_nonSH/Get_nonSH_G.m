function [G_R, G_S, IONC, m0, NN] = Get_nonSH_G(fig,doy,Sites_Info,sate,SDCB_REF,K,M,PG,sate_mark)
%%  estimate satellite and receiver DCBs, ionospheric parameters, et al.
%%  produced from 'Get_MDCB.m' in M_DCB 
% INPUT:
%     fig: group number of products
%     doy: year and doy of year
%     Sites_Info: name and coordinate information of the stations
%     sate: precise coordinates of the satellites
%     SDCB_REF: reference satellite DCBs
%     K,M: order and degree of nonintegral-SH model
%     PG: weight of GPS observations
% OUTPUT:
%     G_R, G_S: estimated receiver and satellite DCBs
%     IONC: ionospheric parameters
%     m0: standard deviation
%     NN: covariance matrices
%% written by Jin R et al., 2012/5/26, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/14
%% --------------------------------------------------------------------------
global sample_num;
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
path_G=['P4/regional/GPS/' doy];
list_gps=dir([path_G '/*.mat']);
G_n_r=length(list_gps);%the number of receivers
%--check the number of each satellite's observations 
gpsnum=sum(sate_mark.gps);
G_PRN=linspace(0,0,gpsnum);
G_S=linspace(0,0,gpsnum);
for i=1:G_n_r
    load([path_G '/' list_gps(i).name],'-mat');
    for j=1:gpsnum
        for k=1:sample_num
            if GPSP4(k,j)~=0
                G_PRN(j)=G_PRN(j)+1;
            end
        end
    end
    clear GPSP4;
end
d_sat=find(G_PRN==0);
if isempty(d_sat)
    G_n_s=gpsnum;
else
    G_n_s=gpsnum-length(d_sat);%the number of satellites
    disp([doy ' PRN ',num2str(d_sat) ,' have no observations.']);
    for k=length(d_sat):-1:1
        gpsx(:,d_sat(k))=[];gpsy(:,d_sat(k))=[];gpsz(:,d_sat(k))=[];
    end
end

%--chose the order of spheric harmonic function
%--LS estimate
num=((K+1)^2-(K-M)*(K-M+1))*fig+G_n_s+G_n_r;
N=zeros(num,num);
U=zeros(num,1);
L=0; sizel=0;
C=linspace(0,0,num);
C(G_n_r+1:G_n_r+G_n_s)=ones(1,G_n_s); 
if G_n_s==gpsnum
    Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy); 
    Wx=-sum(SDCB_REF.gps(index,d_sat));
    G_S(d_sat)=SDCB_REF.gps(index,d_sat);
end
for i=1:G_n_r
    load([path_G '/' list_gps(i).name],'-mat');
    if ~isempty(d_sat)
        for k=length(d_sat):-1:1
            GPSP4(:,d_sat(k))=[];
        end
    end
    site=list_gps(i).name(1:4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    [sN,sl]=Get_Matrix(fig,GPSP4,gpsx,gpsy,gpsz,sx,sy,sz,G_n_r,G_n_s,i,K,M);
    N=N+sN'*sN*PG;
    U=U+sN'*sl*PG;
    % --RMS
    sizel=sizel+length(sl);
    L=L+sl'*sl*PG;
    %-----------
    clear GPSP4;
    disp(['------------>  [',num2str(i),'/',num2str(G_n_r),']  ',num2str(i/G_n_r*100),' % GPS data has been solved!!']);
end

N=N+C'*C;
U=U+C'*Wx;
L=L+Wx'*Wx;
R=pinv(N)*U;
G_R=R(1:G_n_r)*10^9/299792458;
temp=linspace(1,gpsnum,gpsnum);
temp(d_sat)=[];
G_S(temp)=R(G_n_r+1:G_n_r+G_n_s)*10^9/299792458;
IONC=R(G_n_r+G_n_s+1:end);
%--RMS
V=L-R'*U;
f=sizel-num;
m0=sqrt(V/f);
NN=N(G_n_s+G_n_r+1:end,G_n_s+G_n_r+1:end);
%------
end

%% ------------------------get matrix---------------------------
function [MC,l]=Get_Matrix(fig,P4,x,y,z,sx,sy,sz,n_r,n_s,ith,K,M)
MC=[];
l=[];
global sample_num;
num=(K+1)^2-(K-M)*(K-M+1);
[sb,sl]=XYZtoBLH(sx,sy,sz);
figt=sample_num/fig;
for i=1:fig
    for j=1:n_s                %-----------------------j is satellite number
        parfor k=figt*i-(figt-1):figt*i %-------------------------k is epoch number
            if P4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,num*fig+n_s+n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j)*1000,y(k,j)*1000,z(k,j)*1000);
            IPPz=asin(6371000*sin(pi/2-E)/(6371000+450000));%-----SLM 
            t_r=30*(k-1)*pi/43200;
            [b,s]=Get_IPP(E,A,sb,sl,IPPz,t_r);
            M_col(ith)=(-9.52437)*cos(IPPz);   %-----station dcb coefficient
            M_col(n_r+j)=(-9.52437)*cos(IPPz); %---satallite dcb coefficient
            st=num*(i-1)+n_r+n_s+1;
            ed=num*i+n_r+n_s;
            M_col(st:ed)=Get_nonSH(b,s,K,M);
            M_scol=sparse(M_col);
            MC=[MC;M_scol]; 
            l=[l;P4(k,j)*(-9.52437)*cos(IPPz)];
        end
    end
end
end

%% 
function cof_Cap=Get_nonSH(b,s,K,M)

cof_Cap=linspace(0,0,(K+1)^2-(K-M)*(K-M+1));
ms=linspace(s,M*s,M);
i=1;
MM=M;
x=cos(b);
for k=0:K
        P=legendre(k,x);
    if k<M
        M=k;
    end
    for m=0:M
        if m==0
            cof_Cap(i)=P(m+1)*norm(k,m);                    %------------an0
        else
            cof_Cap(i)=P(m+1)*norm(k,m)*cos(ms(m));         %------------anm
            i=i+1;
            cof_Cap(i)=P(m+1)*norm(k,m)*sin(ms(m));         %------------bnm
        end
        i=i+1;
    end
    M=MM;
end
end

%% ------------------------------sub_function--------------------------------
function N=norm(n,m)
if m==0
    N=sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
else
    N=sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
end
end

