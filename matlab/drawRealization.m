function [tauHat, tauTilde, liarMask] = drawRealization(scenario)
%DRAWREALIZATION Sample one round of reports (paper Eq. 14) and surveillance.
%   liarMask(i) = true  -> vehicle i sent a malicious report this round
%   liarMask(i) = false -> vehicle i reported truthfully

N = scenario.N;

u = rand(N,1);
liarMask = u > scenario.lambdaTrue(:);   % truthful w.p. lambdaTrue_i

noise = (2*rand(N,1)-1) .* scenario.eps(:);  % Uniform(-eps_i, eps_i)
tauTilde = scenario.tau(:) + noise;

tauHat = scenario.tau(:);
tauHat(liarMask) = tauHat(liarMask) + scenario.malDelta(liarMask);

end
