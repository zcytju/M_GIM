function Write_ionex(doy,fig,G_R,C_R,E_R,EX_R,R_R,G_S,C_S,E_S,R_S,Pname,VTEC,sate_mark,total_station,list_P4)
%% write ionospheric results to standard IONEX files
% INPUT:
%     doy: year and doy of year
%     fig: group number of products
%     G_R,C_R,E_R,EX_R,R_R: receiver DCBs of different systems
%     G_S,C_S,E_S,R_S: satellite DCBs of different systems
%     Pname: specified name
%     VTEC: ionospheric VTEC results
%     sate_mark: satellite status identification
%     total_station: total number of stations
%     list_P4: P4 observation list for different systems
%% written by Zhou C. et al., 2020/7/15
%% -----------------------------------------------------------------------
figt=24/fig;
interval=figt*3600;
doystr=num2str(doy);
doya=doystr(3:5);doyb=doystr(1:2);
fileID=fopen(['M_Result/',Pname,doya,'0.',doyb,'i'],'w');

sdcbsize=size(G_S,2)+size(C_S,2)+size(E_S,2)+size(R_S,2);
rdcbsize=total_station;
fprintf(fileID,'     1.0            IONOSPHERE MAPS     MIX                 IONEX VERSION / TYPE\n');
dt = datetime(2000+str2num(doyb),1,str2num(doya));
timey=year(dt); timem=month(dt); timed=day(dt);
time_end=(fig-1)*figt;time_h=fix(time_end);time_m=mod(time_end,1)*60;
fprintf(fileID,'%6.0f    %2.0f    %2.0f     0     0     0                        EPOCH OF FIRST MAP  \n',timey,timem,timed);
fprintf(fileID,'%6.0f    %2.0f    %2.0f    %2.0f    %2.0f     0                        EPOCH OF LAST MAP   \n',timey,timem,timed,time_h,time_m);
fprintf(fileID,['%6d                                                      INTERVAL            \n',...
                       '%6d                                                      # OF MAPS IN FILE   \n',...
                       '  COSZ                                                      MAPPING FUNCTION    \n',...
                       '     0.0                                                    ELEVATION CUTOFF    \n',...
                       'combined TEC calculated as weighted mean of input TEC valuesOBSERVABLES USED    \n'],interval,fig);
fprintf(fileID,'%6d                                                      # OF STATIONS       \n',rdcbsize);
fprintf(fileID,'%6d                                                      # OF SATELLITES     \n',sdcbsize);
fprintf(fileID,['  6371.0                                                    BASE RADIUS         \n',...
                       '     2                                                      MAP DIMENSION       \n',...
                       '   450.0 450.0   0.0                                        HGT1 / HGT2 / DHGT  \n']);
latstart=VTEC(1,2); lonstart=VTEC(1,3);
latend=VTEC(end,2); lonend=VTEC(end,3);
fprintf(fileID,'%8.1f%6.1f  -2.5                                        LAT1 / LAT2 / DLAT  \n',latstart,latend);
fprintf(fileID,'%8.1f%6.1f   5.0                                        LON1 / LON2 / DLON  \n',lonstart,lonend);
fprintf(fileID,['    -1                                                      EXPONENT            \n',...
                       'TEC values in 0.1 tec units;  9999, if no value available   COMMENT             \n',...
                       'DCB values in nanoseconds, reference is Sum_of_SatDCBs = 0  COMMENT             \n',...
                       'DIFFERENTIAL CODE BIASES                                    START OF AUX DATA   \n']);
% add satellites DCB
k=1;
for i=1:size(sate_mark.gps,2)
    if sate_mark.gps(i)==1
       fprintf(fileID,'   G%2.0f %9.3f     0.000                                  PRN / BIAS / RMS    \n',i,G_S(k));
       k=k+1;
    else
        continue;
    end
end

k=1;
for i=1:size(sate_mark.bds,2)
    if sate_mark.bds(i)==1
       fprintf(fileID,'   C%2.0f %9.3f     0.000                                  PRN / BIAS / RMS    \n',i,C_S(k));
       k=k+1;
    else
        continue;
    end
end

k=1;
for i=1:size(sate_mark.gal,2)
    if sate_mark.gal(i)==1
       fprintf(fileID,'   E%2.0f %9.3f     0.000                                  PRN / BIAS / RMS    \n',i,E_S(k));
       k=k+1;
    else
        continue;
    end
end

k=1;
for i=1:size(sate_mark.glo,2)
    if sate_mark.glo(i)==1
       fprintf(fileID,'   R%2.0f %9.3f     0.000                                  PRN / BIAS / RMS    \n',i,R_S(k));
       k=k+1;
    else
        continue;
    end
end
% add receiver DCB
list_obs=list_P4.gps;
for i=1:size(G_R,1)
        site=string(list_obs(i).name(1:4));
fprintf(fileID,'   G  %s              %12.3f     0.000              STATION / BIAS / RMS\n',site,G_R(i));
end

list_obs=list_P4.bds;
for i=1:size(C_R,1)
        site=string(list_obs(i).name(1:4));
fprintf(fileID,'   C  %s              %12.3f     0.000              STATION / BIAS / RMS\n',site,C_R(i));
end

list_obs=list_P4.gal;
for i=1:size(E_R,1)
        site=string(list_obs(i).name(1:4));
fprintf(fileID,'   E  %s              %12.3f     0.000              STATION / BIAS / RMS\n',site,E_R(i));
end

list_obs=list_P4.galx;
for i=1:size(EX_R,1)
        site=string(list_obs(i).name(1:4));
fprintf(fileID,'   E  %s              %12.3f     0.000              STATION / BIAS / RMS\n',site,EX_R(i));
end

list_obs=list_P4.glo;
for i=1:size(R_R,1)
        site=string(list_obs(i).name(1:4));
fprintf(fileID,'   R  %s              %12.3f     0.000              STATION / BIAS / RMS\n',site,R_R(i));
end

fprintf(fileID,['DIFFERENTIAL CODE BIASES                                    END OF AUX DATA\n',...     
                       '                                                            END OF HEADER       \n']);
figstart=VTEC(1,1);figend=VTEC(end,1);
fig=(figend-figstart+figt)/figt;
linei=size(VTEC,1)/((figend-figstart+figt)/figt);
linej=linei/((latstart-latend)/2.5+1);
% add TEC data
m=1;
for i=figstart:figt:figend
    hou=floor(i);
    minu=(i-hou)*60;
    fprintf(fileID,'    %2d                                                      START OF TEC MAP    \n',m);
    fprintf(fileID,'%6.0f%6.0f%6.0f%6.0f%6.0f     0                        EPOCH OF CURRENT MAP\n',timey,timem,timed,hou,minu);
    n=1;
    for j=latstart:-2.5:latend
        fprintf(fileID,'%8.1f%6.1f%6.1f   5.0 450.0                            LAT/LON1/LON2/DLON/H\n',j,lonstart,lonend);
        p=1;
        for k=lonstart:5:lonend
            pos=linei*(m-1)+linej*(n-1)+p;
            fprintf(fileID,'%5.0f',VTEC(pos,4)*10);
            if mod(p,16)==0 || p==linej
                fprintf(fileID,'\n');
            end
            p=p+1;
        end
        n=n+1;
    end
    fprintf(fileID,'    %2d                                                      END OF TEC MAP      \n',m);
    m=m+1;
end
% add RMS data
m=1;
for i=figstart:figt:figend
    hou=floor(i);
    minu=(i-hou)*60;
    fprintf(fileID,'    %2d                                                      START OF RMS MAP    \n',m);
    fprintf(fileID,'%6.0f%6.0f%6.0f%6.0f%6.0f     0                        EPOCH OF CURRENT MAP\n',timey,timem,timed,hou,minu);
    n=1;
    for j=latstart:-2.5:latend
        fprintf(fileID,'%8.1f%6.1f%6.1f   5.0 450.0                            LAT/LON1/LON2/DLON/H\n',j,lonstart,lonend);
        p=1;
        for k=lonstart:5:lonend
            pos=linei*(m-1)+linej*(n-1)+p;
            fprintf(fileID,'%5.0f',VTEC(pos,5)*10);
            if mod(p,16)==0 || p==linej
                fprintf(fileID,'\n');
            end
            p=p+1;
        end
        n=n+1;
    end
    fprintf(fileID,'    %2d                                                      END OF RMS MAP      \n',m);
    m=m+1;
end
fprintf(fileID,'                                                            END OF FILE         ');
fclose(fileID);
disp('--------> Ionospheric file has been written!');
end