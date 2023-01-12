function [ ODF_Pseudo ] = getODFPseudo( Icos, d, v, RF, F )
%UNTITLED3 Summary of this function goes here
%   Icos is the Icosahedron mesh generated with IcosahedronMesh.
%   d is the eigenvector matrix of the SPD defining the ellipsoid
%   v is the eigenvector matrix of the SPD defining the ellipsoid
%   R is the rotation matrix from the polar decomposition of the
%   deformation gradient F.
%   F is the deformation gradient

%   Differs from the other approach. Here we do Baroca's 1997 approach but
%   use the sqrt of the stretch instead.

P = v*d*inv(v);
Fapp = F;
F = Fapp*P;
B = F*transpose(F);
[X, b] = eig(B);
pp = X*sqrt(b)*inv(X); pp = pp * 3/trace(pp);
A=zeros(size(Icos));
dv=zeros(size(Icos));
for i=1:length(Icos)
    dv(i,:) = pp*Icos(i,:)';
    A(i,:) = dv(i,:)/norm(dv(i,:));
end

ODF_Pseudo = zeros(length(Icos),1);
ODF_Pseudo = MapToUnitSphere(Icos,A,60);
end 