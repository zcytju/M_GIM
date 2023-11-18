function AreaTEC=Get_areaTEC(fig,lat2,lat1,lon1,lon2,Ref_Data)
%% get VTEC data from IAACs' IONEX files
% INPUT:
%     fig: group number of products
%     lat2, lat1: latitude range
%     lon1, lon2: longitude range
%     Ref_Data: reference IONEX data from IAACs
% OUTPUT:
%     AreaTEC: range limited reference data
%% written by Zhou C. et al., 2020/12/20
%% -------------------------------------------------------
numPD=size(Ref_Data,1);
AreaTEC=[];
figt=24/fig;
for k=0:figt:24
    for i=lat1:-2.5:lat2
        for j=lon1:5:lon2
            for l=1:numPD
                if  k==Ref_Data(l,1)&&i==Ref_Data(l,2)&&j==Ref_Data(l,3)
                    AreaTEC=[AreaTEC;Ref_Data(l,:)];
                end
            end
        end
    end
end
end
