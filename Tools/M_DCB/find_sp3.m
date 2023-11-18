function index=find_sp3(list,doy)
%find related sp3
len=length(list);
doys=linspace(0,0,len);
for i=1:len
    doys(i)=str2double(list(i).name(3:7));
end
index=find(doys==str2double(doy),1);
end 
