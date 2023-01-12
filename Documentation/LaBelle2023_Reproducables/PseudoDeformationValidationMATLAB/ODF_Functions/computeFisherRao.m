% Adam Rauff
% MRL - angiogenesis
% 4/11/2019
% This function computes the distance between 2 Orientation Distribution
% Functions (ODFs) using the Fisher-Rao metric.

% Need to be comparable i.e. same number of points, same normalization
% scheme. Whe

% Goh, Alvina, et al. "A Nonparametric Riemannian Framework for Processing
% High Angular Resolution Diffusion Images (HARDI)" IEEE, 2009

% this function return the distance in radiancs between 2 ODFs
function [ distDeg ] = computeFisherRao(ODF1, ODF2)

% The formulation for this assumes that functions are continuous unit
% vector that sum to 1. This the leads the square root transform to yield
% unit vectors over the unit Hilbert sphere. 
p1 = sqrt(ODF1);
p2 = sqrt(ODF2);


% compute tanget vectors
% distance
distDeg = acosd( dot(p1,p2) );
end

