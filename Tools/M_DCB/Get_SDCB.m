function [DCB_R DCB_S IONC] = Get_SDCB(site,doy,jth,Sites_Info,path_sp3,SDCB_REF,order)
%   Detailed explanation goes here
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
list_sp3=dir([path_sp3 '/*.mat']);
list_obs=dir(['S_P4/' site '/*.mat']);
load(['S_P4/' site '/' list_obs(jth).name],'-mat');
DCB_S=linspace(0,0,32);
%--Test PRN
PRN=linspace(0,0,32);
for j=1:32
    for k=1:2880
        if P4(k,j)~=0
            PRN(j)=PRN(j)+1;
        end
    end
end
d_sat=find(PRN==0);
if isempty(d_sat)
    n_s=32;
else
    n_s=32-length(d_sat);%the number of satellites
    display([site  ' ' doy ' PRN ',num2str(d_sat) ,' have no observations.']);
end
C=linspace(0,0,(order+1)^2*12+n_s+1);
C(2:1+n_s)=ones(1,n_s);
if n_s==32
    Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy);
    Wx=-sum(SDCB_REF.value(index,d_sat));
    DCB_S(d_sat)=SDCB_REF.value(index,d_sat);
end

%--receiver's coordinate

indices=doys==str2double(doy);
index=find(strcmpi(site,stations(indices)), 1);
sx=Coor(index,1);
sy=Coor(index,2);
sz=Coor(index,3);

%--satelliats's coordinate
index2=find_sp3(list_sp3,doy);
load([path_sp3 '/' list_sp3(index2).name],'-mat');

%--ls estimate
[B,l]=Get_Matrix(P4,sx,sy,sz,sate,order);
if ~isempty(d_sat)
    B(:,d_sat+1)=[];
end

BB=[B;C];
L=[l;Wx];
R=BB\L;
DCB_R=R(1)*10^9/299792458;
temp=linspace(1,32,32);
temp(d_sat)=[];
DCB_S(temp)=R(2:1+n_s)*10^9/299792458;
IONC=R(2+n_s:end);

end

%------------------------------sub_function--------------------------------
function [B,l]=Get_Matrix(P4,sx,sy,sz,sate,order)
B=[];l=[];
x=sate.x;y=sate.y;z=sate.z;
[sb,sl]=XYZtoBLH(sx,sy,sz); 
for i=1:12
    for j=1:32                %-----------------------j is satellite number
        for k=240*i-239:240*i %-------------------------k is epoch number
            if P4(k,j)==0
                continue;
            end
            B_col=linspace(0,0,(order+1)^2*12+33);
            [E,A]=Get_EA(sx,sy,sz,x(k,j)*1000,y(k,j)*1000,z(k,j)*1000);
            IPPz=asin(6378137*sin(0.9782*(pi/2-E))/(6378137+506700));%-----MSLM
            t_r=30*(k-1)*pi/43200;
            [b,s]=Get_IPP(E,A,sb,sl,IPPz,t_r);
            B_col(1)=(-9.52437)*cos(IPPz);   
            B_col(1+j)=(-9.52437)*cos(IPPz);
            st=(order+1)^2*(i-1)+34;
            ed=(order+1)^2*i+33;
            B_col(st:ed)=Get_coef(b,s,order);
            B_scol=sparse(B_col);
            B=[B;B_scol]; 
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



