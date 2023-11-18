function PaintData=read_ionex(fig,str)
%% read ionex file to struction
% INPUT:
%     fig: group number of products
%     str: reading mark; 'TEC' for the ionospheric map, and 'RMS' for the RMS map
% OUTPUT:
%     PaintData: final reading results
%% written by Zhou C. et al., 2020/7/15
%% -----------------------------------------------------------------------
[FileName,PathName] = uigetfile('Files_ionex/*.*i','Input ionex files'); % select the ionex files
f = fopen(fullfile(PathName,FileName), 'r'); 
    Data = fscanf(f,'%c'); 
fclose(f);
TECstart = strfind(Data, ['START OF ',str,' MAP']); 
TECend = strfind(Data, ['END OF ',str,' MAP']); 

flag = 0;
PaintData=zeros(1,4);
figt=24/fig;
for i = 1:fig  
    Time=figt*(i-1);
    CurrentMap = Data(TECstart(i): TECend(i)); 
    TECflag = strfind(CurrentMap, 'LAT/LON1/LON2/DLON/H'); %
    for j = 1:71  
        LAT=87.5-2.5*(j-1);
        TECdata = str2num(char(strsplit(strtrim(CurrentMap(TECflag(j)+20: TECflag(j)+390)))));
        for k = 1:73  
            LON=-180+5*(k-1);
            flag = flag + 1;
            PaintData(flag,1:4) = [Time LAT LON TECdata(k,1)*0.1];
        end
    end
end
end

