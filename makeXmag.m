%% makeXmag.m
% Function to select x magnitudes from the specified log normal
% distribution

function x = makeXmag(xtemp,glomExpMax,No,Nx,glomActMu,glomActSig,glomCorr)

if length(glomCorr)==1
    S = glomCorr*(ones(No)-eye(No))+eye(No);
else
    S = glomCorr;
end
valTemp = MvLogNRand(glomActMu*ones(1,No),glomActSig*ones(1,No),Nx,S);

x = xtemp.*valTemp;

% This section just curbs the max value of the glomerular activity
while sum(sum(x>glomExpMax))>0
    x(x>glomExpMax)=lognrnd(glomActMu,glomActSig,size(x(x>glomExpMax)));
end