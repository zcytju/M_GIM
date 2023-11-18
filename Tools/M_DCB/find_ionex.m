function index = find_ionex( list,doy )
%FIND_IONEX Summary of this function goes here
%   Detailed explanation goes here
len=length(list);
doys=linspace(0,0,len);
for i=1:len
    doys(i)=str2double(list(i).name(5:7))+1000*str2double(list(i).name(10:11));
end
index=find(doys==doy,1);
end

