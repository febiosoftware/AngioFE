function [ ODF_Known ] = getODFKnown( Icos, P, F )
%UNTITLED3 Summary of this function goes here
%   Icos is the Icosahedron mesh generated with IcosahedronMesh.m
%   P is the SPD defining the ellipsoid
%   F is the deformation gradient

scale = 3/trace(P);

ODF_Known = zeros(length(Icos),1);
    for i=1:length(Icos)
        ODF_Known(i) = scale*norm(F*P*Icos(i,:)');
    end
ODF_Known = ODF_Known/sum(ODF_Known); 
end

