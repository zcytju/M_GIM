function Sites_Info=read_rinex(r_ipath,r_opath)
%Extract P1,P2,L1,L2 and receivers coordinates.
% ipath: storage path of rinex files
% opath: storage path of L1, L2 ,P1,P2
%--------------------------------------------------------------------------
% assume upper case *O for RINEX files
list_obs=dir([r_ipath '/*.*O']);
if isempty(list_obs)
  % check for lower case
  list_obs=dir([r_ipath '/*.*o']);
end
% if no RINEX files are found, exit with error
if isempty(list_obs)
  error(sprintf('No RINEX files ending in *o or *O found in %s',r_ipath));
end
%list_obs=dir([r_ipath '/*.*O']);
len=length(list_obs);
%creat the folder which shores the L1, L2, P1, P2 observations
if exist(r_opath,'dir')==0 
    mkdir(r_opath);
end  
%--Get P1,P2,L1,L2 and receivers coordinates
Sites_Info.name=cell(1,len);
Sites_Info.doy=linspace(0,0,len)';
Sites_Info.coor=zeros(len,3);
for i=1:len
    obsn=list_obs(i).name;
    doy=obsn(5:7);
    Sites_Info.name{1,i}=obsn(1:4);
    Sites_Info.doy(i)=1000*str2double(obsn(10:11))+str2double(doy);%yyddd
    fname=[r_ipath '/' obsn];
    [obs,coor]=r_rnx(fname); 
    Sites_Info.coor(i,:)=coor;
    %get rid of rinex files without all four types observation
    if all(all(obs.P1==0))||all(all(obs.P2==0))||all(all(obs.L1==0))||all(all(obs.L2==0))
        continue;
    end 
    %save P1,P2,L1,L2 observations.
    name=[obsn(1:4) obsn(10:11) doy '.mat'];
    filename1=[r_opath,'/',name];
    save(filename1,'obs','-mat');
end
save('Sites_Info.mat','Sites_Info','-mat');
end
%------------------------subfunction r_rnx---------------------------------
function [obs,coor]=r_rnx(path)
%-------------------declare variables-----------
obs.P1=zeros(2880,32);
obs.P2=zeros(2880,32);
obs.L1=zeros(2880,32);
obs.L2=zeros(2880,32);     %--------------------------------GPS observables
coor=linspace(0,0,3);
obst=cell(1,11);           %--------------------------GPS observables types 
loc=linspace(0,0,4);       %-------------------target observables' location
obst_n=0;
fid=fopen(path,'r');
while 1
    line=fgetl(fid);
    if ~ischar(line), break, end      %-------the rinex file is null or end
    %-----get receivers coordinates-----------
    if length(line)>78&&strcmp(line(61:79),'APPROX POSITION XYZ')
        coor(1)=str2double(line(1:14));  % X (WGS-84)
        coor(2)=str2double(line(15:28)); % Y (WGS-84)
        coor(3)=str2double(line(29:42)); % Z (WGS-84)
        continue;
    end
    %-----get the GPS observables' types------
    if length(line)>78&&strcmp('# / TYPES OF OBSERV',line(61:79))
        obst_n=str2double(line(5:6));
        if obst_n>9
            for i=1:9
                obst(i)={line(5+6*i:6+6*i)};
            end             
            line=fgetl(fid);
            for i=1:obst_n-9
                obst(i+9)={line(11+6*(i-1):12+6*(i-1))};
            end       
        else
            for i=1:obst_n
               obst(i)={line(11+6*(i-1):12+6*(i-1))};
            end
        end
        for i=1:obst_n
            if strcmp(obst(1,i),'L1')
                loc(1)=i;
                continue;
            end
            if strcmp(obst(1,i),'L2')
                loc(2)=i;
                continue;
            end
            if strcmp(obst(1,i),'P1')
                loc(3)=i;
                continue;
            end
            if strcmp(obst(1,i),'P2')
                loc(4)=i;
                continue;
            end
        end
        if loc(3)==0||loc(4)==0,break,end
        continue; 
    end
    %----start get GPS observables------------
    if(length(line)>72&&strcmp(line(61:73),'END OF HEADER')) 
        while 1                      
            line=fgetl(fid);
            if ~ischar(line)
                break;
            end     %---------------------come to end,break
            %---get epoch number-------- 
            if length(line)>32 &&(strcmp(line(33),'G')||strcmp(line(33),'R')||strcmp(line(33),' '))
                h=str2double(line(11:12));
                m=str2double(line(14:15));
                s=str2double(line(17:26));
                ep=h*120+m*2+s/30+1;   
            else
                continue;
            end 
            if ep~=fix(ep),continue,end
            %---get satellite number----
            nsat=str2double(line(31:32));
            sv_G=linspace(0,0,nsat);
            if nsat>12
                for i=1:12
                    if strcmp(line(30+3*i),'G')
                        sv_G(i)=str2double(line(31+i*3:32+i*3));
                    end
                end
                line=fgetl(fid);
                if nsat<25
                    for i=1:nsat-12
                        if strcmp(line(30+3*i),'G')||strcmp(line(33),' ')
                            sv_G(i+12)=str2double(line(31+i*3:32+i*3));
                        end
                    end
                else
                    for i=1:12
                        if strcmp(line(30+3*i),'G')||strcmp(line(33),' ')
                            sv_G(i+12)=str2double(line(31+i*3:32+i*3));
                        end
                    end
                    line=fgetl(fid);
                    for i=1:nsat-24
                        if strcmp(line(30+3*i),'G')||strcmp(line(33),' ')
                            sv_G(i+12)=str2double(line(31+i*3:32+i*3));
                        end
                    end                  
                end  
            else
                for i=1:nsat
                    if strcmp(line(30+3*i),'G')||strcmp(line(33),' ')
                        sv_G(i)=str2double(line(31+i*3:32+i*3));
                    end
                end
            end
            %---get the observations----  
            for i=1:nsat
                 line=fgetl(fid);
                 obs_temp=linspace(0,0,obst_n);
                 if obst_n>5
                     if sv_G(i)==0
                         line=fgetl(fid);
                         continue;
                     end
                     for j=1:5
                         if length(line)>16*j-3
                             if  strcmp(line(16*j-15:16*j-2),'              ')
                                 continue;
                             end
                              obs_temp(j)=str2double(line(16*j-15:16*j-2));
                         end 
                     end
                     line=fgetl(fid);
                     for j=1:obst_n-5
                         if length(line)>16*j-3
                             if  strcmp(line(16*j-15:16*j-2),'              ')
                                 continue;
                             end
                              obs_temp(j+5)=str2double(line(16*j-15:16*j-2));
                         end
                     end
                 else
                     if sv_G(i)==0,continue,end
                     for j=1:obst_n
                         if length(line)>16*j-3
                             if  strcmp(line(16*j-15:16*j-2),'              ')
                                 continue;
                             end
                             obs_temp(j)=str2double(line(16*j-15:16*j-2));
                         end
                     end
                 end
                 obs.L1(ep,sv_G(i))=obs_temp(loc(1));
                 obs.L2(ep,sv_G(i))=obs_temp(loc(2));
                 obs.P1(ep,sv_G(i))=obs_temp(loc(3));
                 obs.P2(ep,sv_G(i))=obs_temp(loc(4));
            end
        end
    end
end
fclose(fid);
end
%------------------------subfunction Get_R_coor---------------------------------

