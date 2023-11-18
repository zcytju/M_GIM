function VTECresult = Get_VTEC(fig, latlim, lonlim, IONC, N, m0, K, M, lat2, lat1, lon1, lon2)
%% get VTEC value by Spheriacl Harmonic (SH) or nonSH extension.
% INPUT:
%     fig: group number of products
%     latlim, lonlim: latitude and longitude resolution
%     IONC: ionospheric parameters
%     N: covariance matrices
%     m0: standard deviation
%     K: if nargin<8, K denotes the order and degree of SH model;
%         if nargin>8, K denotes the order of nonSH model
%     M: the degree of nonSH model
%     lat2, lat1: latitude range
%     lon1, lon2: longitude range
% OUTPUT:
%     VTECresult: the estimated ionospheric VTEC result
%% written by Zhou C. et al., 2020/7/15
%% -------------------------------------------------------------------
VTECresult=[];
figt=2880/fig;
if nargin<8
    lat1=87.5;lat2=-87.5;
    lon1=-180;lon2=180;
    order=K;
    for i=1:fig
        ep=(i-1)*figt;
        num=(order+1)^2;
        for b=lat1:-latlim:lat2
            for s=lon1:lonlim:lon2
                bb=b/180*pi;
                ss=s/180*pi;
                cof_P=Get_SH(ep,bb,ss,order);
                ionc=IONC((1+num*(i-1)):num*i);
                VTEC=cof_P*ionc;
                inv_N=inv(N(1+num*(i-1):num*i,1+num*(i-1):num*i));
                RMS=m0*sqrt(cof_P*inv_N*cof_P')*10;
                VTECresult=[VTECresult;(i-1)*(24/fig),b,s,VTEC,RMS];
            end
        end
        disp(['--------> [ ',num2str(i),' / ',num2str(fig),' ] global ionospheric data has been calculated!']);
    end
else
    for i=1:fig
        ep=(i-1)*figt;
        num=(K+1)^2-(K-M)*(K-M+1);
        for b=lat1:-latlim:lat2
            for s=lon1:lonlim:lon2
                bb=b/180*pi;
                ss=s/180*pi;
                cof_P=Get_nonSH(ep,bb,ss,K,M);
                ionc=IONC(1+num*(i-1):num*i);
                VTEC=cof_P*ionc;
                inv_N=pinv(N(1+num*(i-1):num*i,1+num*(i-1):num*i));
                RMS=m0*sqrt(cof_P*inv_N*cof_P')*10;
                VTECresult=[VTECresult;(i-1)*(24/fig),b,s,VTEC,RMS];
            end
        end
        disp(['--------> [ ',num2str(i),' / ',num2str(fig),' ] regional ionospheric data has been calculated!']);
    end
end
end
%% 
function cof_capP=Get_nonSH(ep,b,s,K,M)

UT=ep*30/3600;
cof_capP=linspace(0,0,(K+1)^2-(K-M)*(K-M+1));
s=s+(UT-12)*15*pi/180;
ms=linspace(s,M*s,M);
i=1;
MM=M;
x=cos(b);
for k=0:K
        P=legendre(k,x);
    if k<M
        M=k;
    end
    for m=0:M
        if m==0
            cof_capP(i)=P(m+1)*norm(k,m);                    %------------an0
        else
            cof_capP(i)=P(m+1)*norm(k,m)*cos(ms(m));         %------------anm
            i=i+1;
            cof_capP(i)=P(m+1)*norm(k,m)*sin(ms(m));         %------------bnm
        end
        i=i+1;
    end
    M=MM;
end
end

%------------------------------sub_function--------------------------------
function cof_P=Get_SH(ep,b,s,order)
UT=ep*30/3600;
S=s+(UT-12)*15*pi/180;
cof_P=linspace(0,0,(order+1)^2);
ms=linspace(S,order*S,order);
i=1;
x=sin(b);
for n=0:order
    P=legendre(n,x);
    for m=0:n
        if m==0
            cof_P(i)=P(m+1)*norm(n,m);                    %------------an0
        else
            cof_P(i)=P(m+1)*norm(n,m)*cos(ms(m));         %------------anm
            i=i+1;
            cof_P(i)=P(m+1)*norm(n,m)*sin(ms(m));         %------------bnm
        end
        i=i+1;
    end
end
end
%------------------------------sub_function--------------------------------
function N=norm(n,m)
if m==0
    N=sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
else
    N=sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
end
end

