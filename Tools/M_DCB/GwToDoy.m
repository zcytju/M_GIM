function [y,doy]= GwToDoy(GPSweek,day)
%GPSWEEKTODOY Summary of this function goes here
%   Detailed explanation goes here
a=[360,725,1090,1455,1821,2186,2551,2916,3282,3647,4012,4377,...
    4743,5108,5473,5838,6204,6569,6934,7299,7665,8030,8395,...
    8760,9126,9491,9856,10221,10587,10952,11317,11682,12048,...
    12413,12778,13143,13509,13874,14239,14604,14970];
doy=0;
y=1980;
n=GPSweek*7+day;
for i=1:41
    if ~(n>a(1))
        doy=n+6;
        break;
    end
    if ~(n>a(i))
        y=y+i-1;
        doy=n-a(i-1);
        break;
    end
end
end

