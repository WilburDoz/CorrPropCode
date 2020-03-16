%% March 15th - Analytics/Numerics for Correlation Propagation
% In this file we do the numerics to predict the propagation of correlation
% in our model of the olfactory system. The outline is as follows:
%
% i) As a function of class sim (f) find the bulb correlations
% ii) As a function of class sim (f) find the piriform correlations
% iii) Show the off-diagonal covariance as a function of f
% iv) Run simulations to calculate the real propagation of correlation
% v) Plot and compare simulations and theory.

% Lets define all the parameters first
global Sx GlomActSig Sy Nx Ny Sce Sci GlomActMu
Sx = 0.1;
GlomActSig = 0.5;
Sy = 0.062;
Nx = 1000;
Ny = 50000;
Sce = 0.2;
Sci = 0.4;
GlomActMu = 0.1;
%% Find the Bulb Correlations as a function of f
% Here we illustrate the calculated formula for bulb correlation's
% dependence on the class parameter f, and check it makes sense.
fs = [0:0.01:0.99, 0.9955];

for i = 1:length(fs)
    f = fs(i);
    corrbulb(i)= ((f+Sx.*(1-2.*f)).*(f.*(exp(GlomActSig^2)-1)+1)-Sx+f.*Sx^2)/((1-f.*Sx).*(exp(GlomActSig^2)-Sx));
end

figure
subplot(1,3,1)
sgtitle('Some dependencies on f')
plot(fs, corrbulb)
xlabel('Class Similarity Parameter - f')
ylabel('Bulb Correlation - \rho')
xlim([0,1])
ylim([0,1])

%% Next we do the piriform version
corrpiri = zeros(1, length(fs));

% Now start calculating some of the terms that don't change
% First get the deviation of the Neuronal Response
gamma = Nx*Sce*(1+Sce/Sci)*Sx*exp(2*(GlomActMu + GlomActSig^2));

% Find the threshold
theta = 2^(0.5)*gamma *erfcinv(2*Sy);

% Now calculate the mean
FirstMoment = @(x) (1/(gamma*(2*pi)^0.5))*(x-theta).*exp(-x.^2./(2*gamma^2));
Mean = integral(FirstMoment, theta, inf);

% And the second moment
SecondMoment = @(x) (1/(gamma*(2*pi)^0.5))*(x-theta).^2.*exp(-x.^2./(2*gamma^2));
SecMom = integral(SecondMoment, theta, inf);

% Finally we run through a loop to calculate the mixed second moment
for i = 1:length(fs)
    f = fs(i);
    
    F = (f + Sx - 2*f*Sx)/(1-f*Sx)*exp(-GlomActSig^2)*(f*(exp(GlomActSig^2)-1)+1);

    % The function to integrate
    MixedSecMom = @(x,y) (x - theta).*(y-theta).*exp(-(x.^2-x.*y.*2*F+y.^2)./(2*gamma^2*(1-F^2)));
    MixMoment = 1/(2*pi*gamma^2*(1-F^2)^0.5)*integral2(MixedSecMom, theta,inf,theta,inf);
    corrpiri(i) = (MixMoment - Mean^2)/(SecMom - Mean^2);
    
    % Just output probability - to check we're doing it right
    % Should vary between Sy^2 and Sy
    Prob = @(x,y) exp(-(x.^2-x.*y.*2*F+y.^2)./(2*gamma^2*(1-F^2)));
    disp(1/(2*pi*gamma^2*(1-F^2)^0.5)*integral2(Prob, theta, inf,theta,inf))
end

subplot(1,3,2)
plot(fs, corrpiri)
xlabel('Class Similarity Parameter - f')
ylabel('Piriform Correlation - \rho')
xlim([0,1])
ylim([0,1])
%% Plotting the behaviour of big F
Fs = zeros(1,length(fs));

for i = 1:length(fs)
    f = fs(i);
    
    Fs(i) = (f + Sx - 2*f*Sx)/(1-f*Sx)*exp(-GlomActSig^2)*(f*(exp(GlomActSig^2)-1)+1);
end

subplot(1,3,3)
plot(fs, Fs)
xlabel('Class Similarity Parameter - f')
ylabel('Covariance Matrix Element - F')
xlim([0,1])
ylim([0,1])
%% Finally we run the simulations to check with theory
f = [0:0.05:1];
NumExp = length(f);
NumberOdours = 100; % Up this for slower more accurate simulation points

PiriSims = zeros(1, NumExp);
PiriErrors = PiriSims;
OBSims = PiriSims;
OBErrors = PiriSims;
mask = triu(true(NumberOdours),1);

for i = 1:length(f)
    % Choose the f value
    ClassSim = f(i);
    disp(['Class Similarity: ',num2str(ClassSim)])
    
    % Create bulb and Piri
    x = makeOdours(NumberOdours, ClassSim);
    y = makePiriform(x, []);

    % Extract pearson correlations in bulb piriform
    Correl = corr(x);
    Correl = Correl(mask);
    OBSims(i) = mean(Correl);
    OBErrors(i) = std(Correl);
    Correl = corr(y);
    Correl = Correl(mask);
    PiriSims(i) = mean(Correl);
    PiriErrors(i) = std(Correl);
end
%% Finally we plot up the comparison
figure
hold on
errorbar(OBSims, PiriSims, PiriErrors, PiriErrors,OBErrors,OBErrors,'*', 'DisplayName', ['Simulations - Bulb Sparseness: ',num2str(Sx)]);
plot(corrbulb, corrpiri, 'DisplayName','Theory')
xlabel('Bulb Correlation')
ylabel('Piriform Correlation')
title('Propagation of Correlation Through Random Matrix')
legend('Location', 'NorthWest')
xlim([0,1])
ylim([0,1])