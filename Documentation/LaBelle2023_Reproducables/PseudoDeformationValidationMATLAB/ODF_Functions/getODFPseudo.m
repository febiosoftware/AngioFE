%function [ ODF_Pseudo ] = getODFPseudo2( Icos, d, v, RF, F )
function [ ODF_Pseudo ] = getODFPseudo( Icos, P, F )

%calculate the net deformation from a unit sphere
Fn = F*P;
[~, ~, Vn] = poldecomp(Fn);

%normalize values
scale = 3/trace(P);

%calculate the value for each face on the icosahedron
ODF_Pseudo = zeros(length(Icos),1);
    for i=1:length(Icos)
        ODF_Pseudo(i) = scale*norm(Vn*Icos(i,:)');
    end
%normalize CDF to sum to 1
ODF_Pseudo = ODF_Pseudo/sum(ODF_Pseudo); 

end 