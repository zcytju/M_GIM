function VTEC=Get_VTEC_POLY(fig,latlim,lonlim,IONC,m0,N,K,M,lat2,lat1,lon1,lon2,lat0,lon0)
%% get VTEC from polynomial model
% INPUT:
%     fig: group number of products
%     latlim, lonlim: latitude and longitude resolution
%     IONC: ionospheric parameters
%     N: covariance matrices
%     m0: standard deviation
%     K: the order of polynomial model
%     M: the degree of polynomial model
%     lat2, lat1: latitude range
%     lon1, lon2: longitude range
% OUTPUT:
%     VTEC: the estimated ionospheric VTEC result
%% written by Zhou C. et al., 2020/7/15
%% -------------------------------------------------------------------
figt=2880/fig;
VTEC=[];
for i=1:fig
    ep=(i-1)*figt;
    num=(K+1)*(M+1);
    for b=lat1:-latlim:lat2
        for s=lon1:lonlim:lon2
            bb=b/180*pi;
            ss=s/180*pi;
            cof_P=Get_poly(ep,bb,ss,K,M,lat0,lon0);
            ionc=IONC(1+num*(i-1):num*i);
            vtec=cof_P*ionc;
            inv_N=pinv(N(1+num*(i-1):num*i,1+num*(i-1):num*i));
            RMS=m0*sqrt(cof_P*inv_N*cof_P')*10;
            VTEC=[VTEC;(i-1)*(24/fig),b,s,vtec,RMS];
        end
    end
    disp(['--------> [ ',num2str(i),' / ',num2str(fig),' ] regional ionospheric data has been calculated!']);
end
end

%% 
function cof_P=Get_poly(ep,b,s,K,M,lat0,lon0)

UT=ep*30/3600;
cof_P=linspace(0,0,(K+1)*(M+1));
s=s+(UT-12)*15*pi/180;

diff_b=b-lat0;
diff_s=s-lon0;
m=1;
for i=0:K
    for j=0:M
        cof_P(m)=diff_b^i*diff_s^j;
        m=m+1;
    end
end
end
