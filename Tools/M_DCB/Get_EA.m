function [E,A]= Get_EA(sx,sy,sz,x,y,z)
%GET_EL Summary of this function goes here
%   Detailed explanation goes here
[sb,sl]=XYZtoBLH(sx,sy,sz);
T=[-sin(sb)*cos(sl) -sin(sb)*sin(sl) cos(sb);
    -sin(sl)               cos(sl)         0;
    cos(sb)*cos(sl) cos(sb)*sin(sl)  sin(sb)];%transition matrix(XYZ to NEU)
deta_xyz=[x,y,z]-[sx,sy,sz];
NEU=T*(deta_xyz)';
E=atan(NEU(3)/sqrt(NEU(1)*NEU(1)+NEU(2)*NEU(2)));
A=atan(abs(NEU(2)/NEU(1)));
if NEU(1)>0
    if NEU(2)>0
    else
        A=2*pi-A;
    end
else
    if NEU(2)>0
        A=pi-A;
    else
        A=pi+A;
    end 
end
end

