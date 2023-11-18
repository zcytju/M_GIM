function [] = Get_P12(path_obs,path_sp3,Sites_Info,lim,flag)
%--get smoothed P4 observations
%flag=0 : P4 files sorted by doy (estimate based on multi stations)
%flag~=0 : P4 files sorted by sites name ((estimate based on single stations)
%   Detailed explanation goes here
list_obs=dir([path_obs '/*.mat']);
list_sp3=dir([path_sp3 '/*.mat']);
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
if exist('M_P4','dir')==0
    mkdir('M_P4');
end
len=length(list_obs);
for i=1:len
    load([path_obs,'/',list_obs(i).name],'-mat');
    site=list_obs(i).name(1:4);
    doy=list_obs(i).name(5:9);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);   
    if sx==0&sy==0&sz==0
        continue;
    end %get rid of sites without receivers coordinates.
    sp3index=find_sp3(list_sp3,doy);
    load([path_sp3 '/' list_sp3(sp3index).name],'-mat');
    obs=cutobs(sate,sx,sy,sz,obs,lim);
    if all(all(obs.P1==0))||all(all(obs.P2==0))
        continue;
    end
    P4=prepro(obs); 
    if flag==0
        if exist(['M_P4/' doy],'dir')==0
            mkdir(['M_P4/' doy]);
        end
        filenameP4=['M_P4/' doy '/' site doy 'P4.mat'];
        save(filenameP4,'P4','-mat');
    else
        if exist(['S_P4/' site],'dir')==0
            mkdir(['S_P4/' site]);
        end
        filenameP4=['S_P4/' site '/' site doy 'P4.mat'];
        save(filenameP4,'P4','-mat');
    end

    clear P4;
    clear sate;
end
end
%----------------subfunction-----------------------------------------------
function obs=cutobs(sate,sx,sy,sz,obs,lim)
x=sate.x;y=sate.y;z=sate.z;
for i=1:32
    for j=1:2880
        if obs.L1(j,i)==0||obs.L2(j,i)==0||obs.P1(j,i)==0||obs.P2(j,i)==0
            obs.P1(j,i)=0;obs.P2(j,i)=0;obs.L1(j,i)=0;obs.L2(j,i)=0;
            continue;
        end 
        [el,aaa]=Get_EA(sx,sy,sz,x(j,i)*1000,y(j,i)*1000,z(j,i)*1000);
        if el<lim
            obs.P1(j,i)=0;obs.P2(j,i)=0;obs.L1(j,i)=0;obs.L2(j,i)=0;
            continue;
        end
    end
end
end
%----------------subfunction-----------------------------------------------
function P4=prepro(obs)
 %____delete incomplete epoch_____
 P4=zeros(2880,32);
 L4=zeros(2880,32);
 c=299792458;                     %__________________________speed of light
 f1=1575.42*10^6;                 %_________________________________unit:Hz
 f2=1227.6*10^6;                  %_________________________________unit:Hz 
 lamda_w=299792458/(f1-f2);       %____________________wide lane wavelength
 L6=lamda_w*(obs.L1-obs.L2)-(f1*obs.P1+f2*obs.P2)/(f1+f2); %__MW observable 
 Li=obs.L1-f1*obs.L2/f2;
 Nw=L6/lamda_w;                   %_____________________wide lane ambiguity
 for i=1:32  %i is PRN number
     %------divide arc---------------------------
     arc=Get_arc(L6(:,i));
     [arc_n,aaa]=size(arc);
     %----delete arc less than 10 epoches-------
     arc_d=[];
     for j=1:arc_n
         n_epoch=arc(j,2)-arc(j,1);
         if n_epoch<10
             for k=arc(j,1):arc(j,2)
                 obs.P1(k,i)=0;obs.P2(k,i)=0;obs.L1(k,i)=0;obs.L2(k,i)=0;
                 L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
             end
             arc_d=[arc_d,j];
         end
     end
     arc(arc_d,:)=[];
     %----mw detect cycle slip------------------
     [arc_n,aaa]=size(arc);
     %slip=[];
     j=1;
     while j<arc_n+1 %j is arc number 
         %----first epoch check----------
         e=arc(j,1);
         while 1
             if e+1==arc(j,2)||e==arc(j,2)
                 break;
             end
             fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
             firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
             sub=abs(fir-sec);sub2=abs(sec-thi);
             subl=abs(firl-secl);subl2=abs(secl-thil);
             if sub>1||sub2>1||subl>1||subl2>1
                 L6(e,i)=0;obs.L1(e,i)=0;obs.L2(e,i)=0;obs.P1(e,i)=0;obs.P2(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                 e=e+1;
                 arc(j,1)=e;
             else
                 arc(j,1)=e;
                 break;
             end
         end
         %----detect------------------
         if arc(j,2)-arc(j,1)<10
             for k=arc(j,1):arc(j,2)
                 obs.P1(k,i)=0;obs.P2(k,i)=0;obs.L1(k,i)=0;obs.L2(k,i)=0;
                 L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
             end
             arc(j,:)=[];
             arc_n=arc_n-1;
             continue;
         end
         ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
         sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
         sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
         ave_N(1)=Nw(arc(j,1),i);
         sigma2(1)=0;
         sigma(1)=0;
         count=2;
         for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
              ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
              sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
              sigma(count)=sqrt(sigma2(count));
              T=abs(Nw(k+1,i)-ave_N(count));
              I1=abs(Li(k+1,i)-Li(k,i));
              if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                  count=count+1;
                  continue;
              else
                  if k+1==arc(j,2)            %---------------------arc end
                     if k+1-arc(j,1)>10
                         L6(k+1,i)=0;obs.P1(k+1,i)=0;obs.P2(k+1,i)=0;
                         obs.L1(k+1,i)=0;obs.L2(k+1,i)=0;Nw(k+1,i)=0;
                         Li(k,i)=0;arc(j,2)=k;
                     else                     %------delete scatter epoches
                          for l=arc(j,1):k+1
                              L6(l,i)=0;obs.P1(l,i)=0;obs.P2(l,i)=0;
                              obs.L1(l,i)=0;obs.L2(l,i)=0;Nw(l,i)=0;
                              Li(k,i)=0;
                          end
                          arc(j,:)=[];
                          j=j-1;
                          arc_n=arc_n-1;
                     end
                     break;
                  end
                  I2=abs(Li(k+2,i)-Li(k+1,i));
                  if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                      if k+1-arc(j,1)>10
                          arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                          arc_n=arc_n+1;
                      else                    %------delete scatter epoches 
                          for l=arc(j,1):k
                              L6(l,i)=0;obs.P1(l,i)=0;obs.P2(l,i)=0;
                              obs.L1(l,i)=0;obs.L2(l,i)=0;Nw(l,i)=0;
                              Li(k,i)=0;
                          end
                          arc(j,1)=k+1;
                          j=j-1;
                      end
                  else                        %-----------------gross error
                      if k+1-arc(j,1)>10
                          L6(k+1,i)=0;obs.P1(k+1,i)=0;obs.P2(k+1,i)=0;
                          obs.L1(k+1,i)=0;obs.L2(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                          arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                          arc_n=arc_n+1;                     
                      else
                          for l=arc(j,1):k+1  %------delete scatter epoches 
                              L6(l,i)=0;obs.P1(l,i)=0;obs.P2(l,i)=0;
                              obs.L1(l,i)=0;obs.L2(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                          end
                          arc(j,1)=k+2;
                          j=j-1;
                      end
                  end
                  break;
              end  
         end
         j=j+1;
     end
     P4(:,i)=obs.P1(:,i)-obs.P2(:,i);
     L4(:,i)=(c/f1)*obs.L1(:,i)-(c/f2)*obs.L2(:,i);
     %--------smoothing-------------------------
     for j=1:arc_n
         t=2;
         for k=arc(j,1)+1:arc(j,2)
             %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
             P4(k,i)=P4(k,i)/t+(P4(k-1,i)+L4(k-1,i)-L4(k,i))*(t-1)/t;
             t=t+1;
         end
        P4(arc(j,1):arc(j,1)+4,i)=0;
     end
     %--------remove bad P4---------------------
     arc=Get_arc(P4(:,i));
     [arc_n,aaa]=size(arc);
     for j=1:arc_n
         ave=mean(P4(arc(j,1):arc(j,2),i));
         if abs(ave)>10
             for kk=arc(j,1):arc(j,2)
                 P4(kk,i)=0;
             end
         end
     end
 end
end

