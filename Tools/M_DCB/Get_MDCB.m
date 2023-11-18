function [DCB_R DCB_S IONC] = Get_MDCB(doy ,Sites_Info,sate,SDCB_REF,order)
%GET_GDCB Summary of this function goes here
%   Detailed explanation goes here
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
x=sate.x;y=sate.y;z=sate.z;
path_P12=['M_P4/' doy];
list_obs=dir([path_P12 '/*.mat']);
n_r=length(list_obs);%the number of receivers
%--check the number of each satellite's observations 
PRN=linspace(0,0,32);
DCB_S=linspace(0,0,32);
for i=1:n_r
    load([path_P12 '/' list_obs(i).name],'-mat');
    for j=1:32
        for k=1:2880
            if P4(k,j)~=0
                PRN(j)=PRN(j)+1;
            end
        end
    end
    clear P4
end
d_sat=find(PRN==0);
if isempty(d_sat)
    n_s=32;
else
    n_s=32-length(d_sat);%the number of satellites
    display([doy ' PRN ',num2str(d_sat) ,' have no observations.']);
end

%--chose the order of spheric harmonic function
%order=str2double(input('Please input the order of spheric harmonic function (4 order is recommended):','s'));
%--LS estimate
B=[];l=[];
C=linspace(0,0,(order+1)^2*12+n_s+n_r);
C(n_r+1:n_r+n_s)=ones(1,n_s); 
if n_s==32
    Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy); 
    Wx=-sum(SDCB_REF.value(index,d_sat));
    DCB_S(d_sat)=SDCB_REF.value(index,d_sat);
end
for i=1:n_r
    load([path_P12 '/' list_obs(i).name],'-mat');
    site=list_obs(i).name(1:4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    ith=i;
    [sN,sl]=Get_Matrix(P4,x,y,z,sx,sy,sz,n_r,ith,order);
    B=[B;sN];
    l=[l;sl];
    clear P4;
end

if ~isempty(d_sat)
    B(:,d_sat+n_r)=[];
end
BB=[B;C];
L=[l;Wx];
R=BB\L;
DCB_R=R(1:n_r)*10^9/299792458;
temp=linspace(1,32,32);
temp(d_sat)=[];
DCB_S(temp)=R(n_r+1:n_r+n_s)*10^9/299792458;
IONC=R(n_r+n_s+1:end);
end

%------------------------------sub_function--------------------------------
function [M,l]=Get_Matrix(P4,x,y,z,sx,sy,sz,n_r,ith,order)
M=[];
%P=[];
l=[];
[sb,sl]=XYZtoBLH(sx,sy,sz); 
for i=1:12
    for j=1:32                %-----------------------j is satellite number
        for k=240*i-239:240*i %-------------------------k is epoch number
            if P4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,(order+1)^2*12+32+n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j)*1000,y(k,j)*1000,z(k,j)*1000);
            IPPz=asin(6378137*sin(0.9782*(pi/2-E))/(6378137+506700));%-----MSLM
            t_r=30*(k-1)*pi/43200;
            [b,s]=Get_IPP(E,A,sb,sl,IPPz,t_r);
            M_col(ith)=(-9.52437)*cos(IPPz);   %-----station dcb coefficient
            M_col(n_r+j)=(-9.52437)*cos(IPPz); %---satallite dcb coefficient
            st=(order+1)^2*(i-1)+n_r+33;
            ed=(order+1)^2*i+n_r+32;
            M_col(st:ed)=Get_coef(b,s,order);
            M_scol=sparse(M_col);
            M=[M;M_scol]; 
            l=[l;P4(k,j)*(-9.52437)*cos(IPPz)];
        end
    end
end
end
%------------------------------sub_function--------------------------------
function cof_P=Get_coef(b,s,order)
cof_P=linspace(0,0,(order+1)^2);
% ms=[s,2*s,3*s,4*s,5*s,6*s,7*s,8*s,9*s,10*s,11*s,12*s,13*s,14*s,15*s];
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



