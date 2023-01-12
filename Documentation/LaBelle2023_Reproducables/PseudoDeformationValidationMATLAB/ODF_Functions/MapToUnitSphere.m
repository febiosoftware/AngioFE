function [ ODF ] = MapToUnitSphere(Ncart, POIs, q)    
%    minD = zeros(length(Ncart),2);
    ODF = zeros(length(Ncart),1);
    for i = 1:size(Ncart,1)

        pt = Ncart(i,:);

        % distance metric is geodesic on unit sphere - angle
        tempMat = ones(size(Ncart));
        tempMat(:,1) = tempMat(:,1)*pt(1);
        tempMat(:,2) = tempMat(:,2)*pt(2);
        tempMat(:,3) = tempMat(:,3)*pt(3);
        Dist = acos(dot(tempMat,POIs,2));
        Dist = rad2deg(Dist);
        clear tempMat;
        ODF(i) = ODF(i)+length(find(Dist < q));
    end
    ODF = ODF/sum(ODF);
end
