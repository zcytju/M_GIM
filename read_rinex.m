 function Sites_Info=read_rinex(r_ipath,r_opath)
%% extract P1,P2,L1,L2 and receivers coordinates.
%% produced from 'read_rinex.m' in M_DCB
% INPUT:
%     r_ipath: storage path of rinex files
%     r_opath: storage path of L1, L2 ,P1,P2 observation files
% OUTPUT:
%     Sites_Info: name and coordinate information of the stations
%% written by Jin R et al., 2012/5/20, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/12
%% -----------------------------------------------------------------------
% assume upper case *O for RINEX files
list_obs1=dir([r_ipath '/*.rnx']);
list_obs2=dir([r_ipath '/*.*o']);
list_obs=[list_obs1;list_obs2];
% if no RINEX files are found, exit with error
if isempty(list_obs)
  error(sprintf('No RINEX files ending in *o or *rnx found in %s',r_ipath));
end
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
    if length(obsn)>25
        doy=obsn(15:19);
        Sites_Info.name{1,i}=upper(obsn(1:4));
        name=[upper(obsn(1:4)) doy '.mat'];
        if exist([r_opath,'/',doy],'dir')==0
            mkdir([r_opath,'/',doy]);
        end
        filename1=[r_opath,'/',doy,'/',name];
    else
        doy=[obsn(10:11),obsn(5:7)];
        Sites_Info.name{1,i}=upper(obsn(1:4));
        name=[upper(obsn(1:4)) doy '.mat'];
        if exist([r_opath,'/',doy],'dir')==0
            mkdir([r_opath,'/',doy]);
        end
        filename1=[r_opath,'/',doy,'/',name];
    end
    
    Sites_Info.doy(i)=str2double(doy);%yyddd
    fname=[r_ipath '/' obsn];
    if isfile(filename1)
        % data were previously loaded and saved
        load(filename1,'coor');
    else
        [obs,coor]=read_rnx3(fname); 
        %save P1,P2,L1,L2 observations.
        if  ~isempty(obs) %~isnan(coor) &&
            save(filename1,'obs','coor','-mat');
        end
    end
    Sites_Info.coor(i,:)=coor;
    disp(['-------- [ ',num2str(i),' / ',num2str(len),' ] obs files are read !']);
end
% save('Sites_Info.mat','Sites_Info','-mat');
end
%% ------------------------read rinex 3.02-----------------------------
function [obs,coor]=read_rnx3(path)
%% read observations from RINEX 3 files
% INPUT: 
% path: path of observation files
% OUTPUT: 
% obs: struct of observation files
% coor: station coordinates
%% ------------------------------------------------------------
gpsmark=1;glomark=1;
galmark=1;bdsmark=1;
galxmark=1;
fid=fopen(path,'r');
while 1
    line=fgetl(fid);
    if ~ischar(line), break, end      %-------the rinex file is null or end
    if length(line)>79 && strcmp(line(61:80),'RINEX VERSION / TYPE')
        if strcmp(line(6),'3')==0
            obs=[];coor=[0,0,0];
            break;
        end
    end
    %-----get receivers coordinates-----------
    if length(line)>78 && strcmp(line(61:79),'APPROX POSITION XYZ')
        coor(1)=str2double(line(1:14));  % X (WGS-84)
        coor(2)=str2double(line(15:28)); % Y (WGS-84)
        coor(3)=str2double(line(29:42)); % Z (WGS-84)
        continue;
    end
    %% -----get the GPS observables' types------
    if length(line)>78 && strcmp('G',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        mark.gps=1;
        gpsloc=zeros(1,4);
        obsgps_n=str2double(line(4:6));
        if obsgps_n>39
            for i=1:13
                gpst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                gpst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                gpst(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgps_n-39
                gpst(i+39)={line(4+4*i:6+4*i)};
            end
        elseif obsgps_n>26 && obsgps_n<=39
            for i=1:13
                gpst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                gpst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgps_n-26
                gpst(i+26)={line(4+4*i:6+4*i)};
            end
        elseif obsgps_n>13 && obsgps_n<=26
            for i=1:13
                gpst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgps_n-13
                gpst(i+13)={line(4+4*i:6+4*i)};
            end
        else
            for i=1:obsgps_n
                gpst(i)={line(4+4*i:6+4*i)};
            end
        end
        for i=1:obsgps_n
            if (strcmp(gpst(1,i),'L1C')) || (strcmp(gpst(1,i),'L1 ')) %L1C -> L1
                gpsloc(1)=i;
                continue;
            end
            if strcmp(gpst(1,i),'L2W') || strcmp(gpst(1,i),'L2 ') %L2W -> L2
                gpsloc(2)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C1W') %C1W -> P1
                gpsloc(3)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2W') %C2W -> P2
                gpsloc(4)=i;
                continue;
            end
        end
        if gpsloc(3)==0 || gpsloc(4)==0 || gpsloc(1)==0 || gpsloc(2)==0
            gpsmark=0;markf.gps=0;
        end
        continue; 
    end
     %% -----get the GLONASS observables' types------
    if length(line)>78 && strcmp('R',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        mark.glo=1;
        gloloc=zeros(1,4);
        obsglo_n=str2double(line(4:6));
        if obsglo_n>39
            for i=1:13
                glot(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                glot(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                glot(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsglo_n-39
                glot(i+39)={line(4+4*i:6+4*i)};
            end
        elseif obsglo_n>26 && obsglo_n<=39
            for i=1:13
                glot(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                glot(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsglo_n-26
                glot(i+26)={line(4+4*i:6+4*i)};
            end
        elseif obsglo_n>13 && obsglo_n<=26
            for i=1:13
                glot(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsglo_n-13
                glot(i+13)={line(4+4*i:6+4*i)};
            end
        else
            for i=1:obsglo_n
                glot(i)={line(4+4*i:6+4*i)};
            end
        end
        for i=1:obsglo_n
            if strcmp(glot(1,i),'L1P') || strcmp(glot(1,i),'L1 ') %L1P -> L1
                gloloc(1)=i;
                continue;
            end
            if strcmp(glot(1,i),'L2P') || strcmp(glot(1,i),'L2 ') %L2P -> L2
                gloloc(2)=i;
                continue;
            end
            if strcmp(glot(1,i),'C1P') %C1P -> P1
                gloloc(3)=i;
                continue;
            end
            if strcmp(glot(1,i),'C2P') %C2P -> P2
                gloloc(4)=i;
                continue;
            end
        end
        if gloloc(3)==0 || gloloc(4)==0 || gloloc(1)==0 || gloloc(2)==0
            glomark=0;markf.glo=0;
        end
        continue; 
    end
     %% -----get the GALILEO observables' types------
    if length(line)>78 && strcmp('E',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        mark.gal=1;
        obsgal_n=str2double(line(4:6));
        if obsgal_n>39
            for i=1:13
                galt(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                galt(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                galt(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgal_n-39
                galt(i+39)={line(4+4*i:6+4*i)};
            end
        elseif obsgal_n>26 && obsgal_n<=39
            for i=1:13
                galt(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                galt(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgal_n-26
                galt(i+26)={line(4+4*i:6+4*i)};
            end
        elseif obsgal_n>13 && obsgal_n<=26
            for i=1:13
                galt(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgal_n-13
                galt(i+13)={line(4+4*i:6+4*i)};
            end
        else
            for i=1:obsgal_n
                galt(i)={line(4+4*i:6+4*i)};
            end
        end
        for i=1:obsgal_n
            if strcmp(galt(1,i),'L1C') || strcmp(galt(1,i),'L1 ')  % L1C  -> L1
                galloc(1)=i;
            elseif strcmp(galt(1,i),'L1X') % L1X -> L1
                galxloc(1)=i;
                continue;
            end
            if strcmp(galt(1,i),'L5Q') || strcmp(galt(1,i),'L5 ')  %L5Q -> L5
                galloc(2)=i;
            elseif strcmp(galt(1,i),'L5X') % L5X -> L5
                galxloc(2)=i;
                continue;
            end
            if strcmp(galt(1,i),'C1C') || strcmp(galt(1,i),'C1 ')  % C1C  -> C1
                galloc(3)=i;
            elseif strcmp(galt(1,i),'C1X') % C1X -> C1
                galxloc(3)=i;
                continue;
            end
            if strcmp(galt(1,i),'C5Q') || strcmp(galt(1,i),'C5 ')  %C5Q -> C5
                galloc(4)=i;
            elseif strcmp(galt(1,i),'C5X') % C5X -> C5
                galxloc(4)=i;
                continue;
            end
        end
        if exist('galloc','var')
            galloc(1,length(galloc)+1:4)=0;
            if galloc(3)==0 || galloc(4)==0 || galloc(1)==0 || galloc(2)==0
                galmark=0;markf.gal=0;
            end
        elseif exist('galxloc','var')
            galxloc(1,length(galxloc)+1:4)=0;
            if galxloc(3)==0 || galxloc(4)==0 || galxloc(1)==0 || galxloc(2)==0
                galxmark=0;markf.galx=0;
            end
        end
        continue; 
    end    
     %% -----get the BDS observables' types------
    if length(line)>78 && strcmp('C',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        mark.bds=1;
        bdsloc=zeros(1,4);
        obsbds_n=str2double(line(4:6));
        if obsbds_n>39
            for i=1:13
                bdst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                bdst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                bdst(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-39
                bdst(i+39)={line(4+4*i:6+4*i)};
            end
        elseif obsbds_n>26 && obsbds_n<=39
            for i=1:13
                bdst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                bdst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-26
                bdst(i+26)={line(4+4*i:6+4*i)};
            end
        elseif obsbds_n>13 && obsbds_n<=26
            for i=1:13
                bdst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-13
                bdst(i+13)={line(4+4*i:6+4*i)};
            end
        else
            for i=1:obsbds_n
                bdst(i)={line(4+4*i:6+4*i)};
            end
        end
        for i=1:obsbds_n
            if strcmp(bdst(1,i),'L2I') %L2I -> L2
                bdsloc(1)=i;
                continue;
            end
            if strcmp(bdst(1,i),'L7I') %L7I -> L7
                bdsloc(2)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C2I') %C2I -> C2
                bdsloc(3)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C7I') %C7I -> C7
                bdsloc(4)=i;
                continue;
            end
        end
        if bdsloc(3)==0 || bdsloc(4)==0 || bdsloc(1)==0 || bdsloc(2)==0
            bdsmark=0;markf.bds=0;
        end
        continue; 
    end    
    
    %% 
    if exist('mark','var') && exist('markf','var')
        field=fieldnames(mark);
        fieldf=fieldnames(markf);
        if length(field)==length(fieldf)
            obs=[];break;
        end
    end
        
    %% ----start getting observables------------
    if(length(line)>72&&strcmp(line(61:73),'END OF HEADER')) 
        while 1                      
            line=fgetl(fid);
            if ~ischar(line)
                break;
            end     %---------------------come to end,break
            %---get epoch number-------- 
            if length(line)>32 && strcmp(line(1),'>')
                h=str2double(line(14:15));
                m=str2double(line(17:18));
                s=str2double(line(20:29));
                ep=h*120+m*2+s/30+1;   
            end 
            if ep~=fix(ep),continue,end

            if gpsmark~=0 && exist('gpsloc','var')
                if length(line)>max(gpsloc)*16 && strcmp(line(1),'G')
                    sv_G=str2double(line(2:3));
                    obs.GPSL1C(ep,sv_G)=str2double(line(16*gpsloc(1)-12:16*gpsloc(1)+1));
                    obs.GPSL2W(ep,sv_G)=str2double(line(16*gpsloc(2)-12:16*gpsloc(2)+1));
                    obs.GPSC1W(ep,sv_G)=str2double(line(16*gpsloc(3)-12:16*gpsloc(3)+1));
                    obs.GPSC2W(ep,sv_G)=str2double(line(16*gpsloc(4)-12:16*gpsloc(4)+1));
                end
            end
            if glomark~=0 && exist('gloloc','var')
                if length(line)>max(gloloc)*16 && strcmp(line(1),'R')
                    sv_R=str2double(line(2:3));
                    obs.GLOL1P(ep,sv_R)=str2double(line(16*gloloc(1)-12:16*gloloc(1)+1));
                    obs.GLOL2P(ep,sv_R)=str2double(line(16*gloloc(2)-12:16*gloloc(2)+1));
                    obs.GLOC1P(ep,sv_R)=str2double(line(16*gloloc(3)-12:16*gloloc(3)+1));
                    obs.GLOC2P(ep,sv_R)=str2double(line(16*gloloc(4)-12:16*gloloc(4)+1));
                end
            end
            if galmark~=0 && exist('galloc','var')
                if length(line)>max(galloc)*16 && strcmp(line(1),'E')
                    sv_E=str2double(line(2:3));
                    obs.GALL1C(ep,sv_E)=str2double(line(16*galloc(1)-12:16*galloc(1)+1));
                    obs.GALL5Q(ep,sv_E)=str2double(line(16*galloc(2)-12:16*galloc(2)+1));
                    obs.GALC1C(ep,sv_E)=str2double(line(16*galloc(3)-12:16*galloc(3)+1));
                    obs.GALC5Q(ep,sv_E)=str2double(line(16*galloc(4)-12:16*galloc(4)+1));
                end
            end
            if galxmark~=0 && exist('galxloc','var')
                if length(line)>max(galxloc)*16 && strcmp(line(1),'E')
                    sv_E=str2double(line(2:3));
                    obs.GALL1X(ep,sv_E)=str2double(line(16*galxloc(1)-12:16*galxloc(1)+1));
                    obs.GALL5X(ep,sv_E)=str2double(line(16*galxloc(2)-12:16*galxloc(2)+1));
                    obs.GALC1X(ep,sv_E)=str2double(line(16*galxloc(3)-12:16*galxloc(3)+1));
                    obs.GALC5X(ep,sv_E)=str2double(line(16*galxloc(4)-12:16*galxloc(4)+1));
                end
            end
            if bdsmark~=0 && exist('bdsloc','var')
                if length(line)>max(bdsloc)*16 && strcmp(line(1),'C')
                    sv_C=str2double(line(2:3));
                    obs.BDSL2I(ep,sv_C)=str2double(line(16*bdsloc(1)-12:16*bdsloc(1)+1));
                    obs.BDSL7I(ep,sv_C)=str2double(line(16*bdsloc(2)-12:16*bdsloc(2)+1));
                    obs.BDSC2I(ep,sv_C)=str2double(line(16*bdsloc(3)-12:16*bdsloc(3)+1));
                    obs.BDSC7I(ep,sv_C)=str2double(line(16*bdsloc(4)-12:16*bdsloc(4)+1));
                end
            end
        end
    end
end
fclose(fid);
end