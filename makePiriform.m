%% makePiriform.m
% Program to create a cortical representation by passing an earlier rep
% through the a random matrix.
%
% INPUT:
% y  - input representation
% th - threshold
% J  - Transition matrix
%
% OUTPUT:
% z - output representation
% J - Transistion matrix, to be used again

%%
function [y, J] = makePiriform(x, J)

global Ny Sce Sci Nx Sy

NumberOdours = size(x, 2);

makeJ = false;
if ~exist('G','var') || isempty(J)
    makeJ = true;
end

if makeJ
    % Now we start making the random J
    J = rand(Ny, Nx);
    J(J<Sci) = -Sce/Sci;
    J(J>(1-Sce)) = 1;
    J(J>0 & J~=1) = 0;
end

% Propagate through
y = J*x;

temp = sort(reshape(y,Ny*NumberOdours,1),'descend');
loc = round(Ny*Sy*NumberOdours); % Location of the last active neuron we want to choose active
theta = temp(loc);
y = y-theta;
y(y<0)=0;
end