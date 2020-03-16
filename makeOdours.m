%% makeOdours.m
% Program to create a glomerular representation of a set of random odours
% Inputs:   No - the number of of odours to be simulated
%           f - No by No matrix decribing odour correlation with one and
%           other.
%           e.g
%           [1 0.7 0
%           0.7 1  0        = 2 highly correlated odours with one non odour
%           0   0  1]
%           If a single number is passes rather than a matrix it is assumed
%           to be modelling a class with corr f.
% Outputs:  x - matrix of glomerular response to each No odour (Nx by No)
%           similarity - vector of similarity to referrence odour which is
%           first - not particularly useful.

%% Calculate the glomerular representation of a set of, potentially correlated, odours.
function [x,similarity] = makeOdours(No, f)

global Sx Nx glomActMu glomActSig

% Here we create the setup if we are only simulating one class therefore No
% by No rule is broken, instead all odours have f overlap.
if length(f) == 1
    cVal = (1-f)*Sx/(1 - f*Sx);
    pMat = [ones( round(f*Sx*Nx), No); cVal*ones( Nx-round(f*Sx*Nx), No)];
    similarity = f*ones(1,No);
    
else
    pMat = nan(Nx,No);
    similarity = nan(1,No);
    zeroclassmarker = 0;
    for j=1:No
        if j==1
            Classf = f(j+1,j);
        elseif j==No
            Classf = f(j-1, j);
        else
            Classf = f(j-1,j);
            if Classf == 0 && ~zeroclassmarker
                Classf2 = f(j, j+1);
                if Classf2 == 0
                   zeroclassmarker = 1;
                else
                   Classf = Classf2;
                end
                
            elseif zeroclassmarker
                Classf2 = f(j, j+1);
                if Classf2
                    zeroclassmarker = 0;
                    Classf = Classf2;
                end
            end
        end
        cVal = (1-Classf)*Sx/(1 - Classf*Sx);
        pMat(:,j) = [ones( round(Classf*Sx*Nx), 1); cVal*ones( Nx-round(Classf*Sx*Nx), 1)];
        similarity(j) = Classf;
    end
end
xtmp = rand(Nx,No) < pMat;
x = makeXmag(xtmp,inf,No,Nx,glomActMu,glomActSig,f);
end

