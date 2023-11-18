function [B,L] = XYZtoBLH(X,Y,Z)
%WGS84坐标转换到大地经纬度
%   Detailed explanation goes here
a=6378137;
e2=0.0066943799013;
L=atan(abs(Y/X));
if Y>0
    if X>0
    else
        L=pi-L;
    end
else
    if X>0
        L=2*pi-L;
    else
        L=pi+L;
    end
end
B0=atan(Z/sqrt(X^2+Y^2));
while 1
    N=a/sqrt(1-e2*sin(B0)*sin(B0));
    H=Z/sin(B0)-N*(1-e2);
    B=atan(Z*(N+H)/(sqrt(X^2+Y^2)*(N*(1-e2)+H)));
    if abs(B-B0)<1e-6;break;end
    B0=B;
end
N=a/sqrt(1-e2*sin(B)*sin(B));
end

