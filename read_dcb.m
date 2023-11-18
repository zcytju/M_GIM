function [SDCB_REF,Sites_Info]=read_dcb( i_ipath,Sites_Info,sate_mark)
%% get the reference DCB values, produced from 'read_ionex.m' in M_DCB
% INPUT:
%     i_ipath: storage path of *.BSX files
%     Sites_Info: name and coordinate information of the stations
%     sate_mark: satellite status identification
% OUTPUT:
%     SDCB_REF: reference satellite DCBs
%     Sites_Info: updated name and coordinate information of the stations
%% written by Jin R et al., 2012/4/25, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/12
%% ---------------------------------------------------------------------------
doys=unique(Sites_Info.doy);
n_d=length(doys);
len=length(Sites_Info.name);
Sites_Info.GPS_RDCB_REF=linspace(0,0,len);
Sites_Info.GLO_RDCB_REF=linspace(0,0,len);
Sites_Info.BDS_RDCB_REF=linspace(0,0,len);
Sites_Info.GAL_RDCB_REF=linspace(0,0,len);
Sites_Info.GALX_RDCB_REF=linspace(0,0,len);
SDCB_REF.gps=zeros(n_d,32);
SDCB_REF.glo=zeros(n_d,24);
SDCB_REF.bds=zeros(n_d,16);
SDCB_REF.gal=zeros(n_d,36);
SDCB_REF.galx=zeros(n_d,36);
SDCB_REF.doy=linspace(0,0,n_d);

for i=1:n_d
    sdoy=doys(i);
    index2= Sites_Info.doy==doys(i);
    [GPS_DCB_rec,SDCB_REF.gps(i,:)]=r_gps_dcb([i_ipath '/CAS0MGXRAP_20' num2str(sdoy) '0000_01D_01D_DCB.BSX'],Sites_Info.name(index2));
    [GLO_DCB_rec,SDCB_REF.glo(i,:)]=r_glo_dcb([i_ipath '/CAS0MGXRAP_20' num2str(sdoy) '0000_01D_01D_DCB.BSX'],Sites_Info.name(index2));
    [BDS_DCB_rec,SDCB_REF.bds(i,:)]=r_bds_dcb([i_ipath '/CAS0MGXRAP_20' num2str(sdoy) '0000_01D_01D_DCB.BSX'],Sites_Info.name(index2));
    [GAL_DCB_rec,SDCB_REF.gal(i,:)]=r_gal_dcb([i_ipath '/CAS0MGXRAP_20' num2str(sdoy) '0000_01D_01D_DCB.BSX'],Sites_Info.name(index2));
    [GALX_DCB_rec,SDCB_REF.galx(i,:)]=r_galx_dcb([i_ipath '/CAS0MGXRAP_20' num2str(sdoy) '0000_01D_01D_DCB.BSX'],Sites_Info.name(index2));
    Sites_Info.GPS_RDCB_REF(index2)=GPS_DCB_rec;
    Sites_Info.GLO_RDCB_REF(index2)=GLO_DCB_rec;
    Sites_Info.BDS_RDCB_REF(index2)=BDS_DCB_rec;
    Sites_Info.GAL_RDCB_REF(index2)=GAL_DCB_rec;
    Sites_Info.GALX_RDCB_REF(index2)=GALX_DCB_rec;
        
    SDCB_REF.gps=SDCB_REF.gps(1:size(sate_mark.gps,2));
    gpsdelete=find(sate_mark.gps==0);
    if ~isnan(gpsdelete)
        k=length(gpsdelete);
        for j=k:-1:1
            SDCB_REF.gps(:,gpsdelete(j))=[];
        end
    end
    
    SDCB_REF.bds=SDCB_REF.bds(1:size(sate_mark.bds,2));
    bdsdelete=find(sate_mark.bds==0);
    if ~isnan(bdsdelete)
        k=length(bdsdelete);
        for j=k:-1:1
            SDCB_REF.bds(:,bdsdelete(j))=[];
        end
    end
    
    SDCB_REF.glo=SDCB_REF.glo(1:size(sate_mark.glo,2));
    glodelete=find(sate_mark.glo==0);
    if ~isnan(glodelete)
        k=length(glodelete);
        for j=k:-1:1
            SDCB_REF.glo(:,glodelete(j))=[];
        end
    end
    
    SDCB_REF.gal=SDCB_REF.gal(1:size(sate_mark.gal,2));
    galdelete=find(sate_mark.gal==0);
    if ~isnan(galdelete)
        k=length(galdelete);
        for j=k:-1:1
            SDCB_REF.gal(:,galdelete(j))=[];
            SDCB_REF.galx(:,galdelete(j))=[];
        end
    end
    
    SDCB_REF.doy(i)=doys(i);
end
end

%% ------------------------------sub function-------------------------------
function [DCB_rec,DCB_sat]=r_galx_dcb(fpath,sites)
fid=fopen(fpath,'r');
n_r=length(sites);
DCB_rec=linspace(0,0,n_r);

while 1
    line=fgetl(fid);
    if ~ischar(line), break, end
    %--satellites' DCB
    if length(line)>100 && strcmpi(line(1:7),' DSB  E') && strcmpi(line(26:33),'C1X  C5X') && strcmpi(line(16:19),'    ')
        prn=str2double(line(13:14));
        DCB_sat(prn)=str2double(line(83:91));
        continue;
    end
    %--receivers' DCB
    if length(line)>100 && strcmpi(line(1:12),' DSB  E    E') && strcmpi(line(26:33),'C1X  C5X')
        index=find(strcmpi(line(16:19),sites), 1);
        if ~isempty(index)
            DCB_rec(index)=str2double(line(83:91));
        end
        continue;
    end
end
fclose(fid);
end
%% ------------------------------sub function-------------------------------
function [DCB_rec,DCB_sat]=r_gal_dcb(fpath,sites)
fid=fopen(fpath,'r');
n_r=length(sites);
DCB_rec=linspace(0,0,n_r);

while 1
    line=fgetl(fid);
    if ~ischar(line), break, end
    %--satellites' DCB
    if length(line)>100 && strcmpi(line(1:7),' DSB  E') && strcmpi(line(26:33),'C1C  C5Q') && strcmpi(line(16:19),'    ')
        prn=str2double(line(13:14));
        DCB_sat(prn)=str2double(line(83:91));
        continue;
    end
    %--receivers' DCB
    if length(line)>100 && strcmpi(line(1:12),' DSB  E    E') && strcmpi(line(26:33),'C1C  C5Q')
        index=find(strcmpi(line(16:19),sites), 1);
        if ~isempty(index)
            DCB_rec(index)=str2double(line(83:91));
        end
        continue;
    end
end
fclose(fid);
end
%% ------------------------------sub function-------------------------------
function [DCB_rec,DCB_sat]=r_bds_dcb(fpath,sites)
fid=fopen(fpath,'r');
n_r=length(sites);
DCB_rec=linspace(0,0,n_r);

while 1
    line=fgetl(fid);
    if ~ischar(line), break, end
    %--satellites' DCB
    if length(line)>100 && strcmpi(line(1:7),' DSB  C') && strcmpi(line(26:33),'C2I  C7I') && strcmpi(line(16:19),'    ')
        prn=str2double(line(13:14));
        DCB_sat(prn)=str2double(line(83:91));
        continue;
    end
    %--receivers' DCB
    if length(line)>100 && strcmpi(line(1:12),' DSB  C    C') && strcmpi(line(26:33),'C2I  C7I')
        index=find(strcmpi(line(16:19),sites), 1);
        if ~isempty(index)
            DCB_rec(index)=str2double(line(83:91));
        end
        continue;
    end
end
fclose(fid);
end
%% ------------------------------sub function-------------------------------
function [DCB_rec,DCB_sat]=r_glo_dcb(fpath,sites)
fid=fopen(fpath,'r');
n_r=length(sites);
DCB_rec=linspace(0,0,n_r);

while 1
    line=fgetl(fid);
    if ~ischar(line), break, end
    %--satellites' DCB
    if length(line)>100 && strcmpi(line(1:7),' DSB  R') && strcmpi(line(26:33),'C1P  C2P') && strcmpi(line(16:19),'    ')
        prn=str2double(line(13:14));
        DCB_sat(prn)=str2double(line(83:91));
        continue;
    end
    %--receivers' DCB
    if length(line)>100 && strcmpi(line(1:12),' DSB  R    R') && strcmpi(line(26:33),'C1P  C2P')
        index=find(strcmpi(line(16:19),sites), 1);
        if ~isempty(index)
            DCB_rec(index)=str2double(line(83:91));
        end
        continue;
    end
end
fclose(fid);
end
%% ------------------------------sub function-------------------------------
function [DCB_rec,DCB_sat]=r_gps_dcb(fpath,sites)
fid=fopen(fpath,'r');
n_r=length(sites);
DCB_rec=linspace(0,0,n_r);
DCB_sat=linspace(0,0,32);

while 1
    line=fgetl(fid);
    if ~ischar(line), break, end
    %--satellites' DCB
    if length(line)>100 && strcmpi(line(1:7),' DSB  G') && strcmpi(line(26:33),'C1W  C2W') && strcmpi(line(16:19),'    ')
        prn=str2double(line(13:14));
        DCB_sat(prn)=str2double(line(83:91));
        continue;
    end
    %--receivers' DCB
    if length(line)>100 && strcmpi(line(1:12),' DSB  G    G') && strcmpi(line(26:33),'C1W  C2W')
        index=find(strcmpi(line(16:19),sites), 1);
        if ~isempty(index)
            DCB_rec(index)=str2double(line(83:91));
        end
        continue;
    end
end
fclose(fid);
end
