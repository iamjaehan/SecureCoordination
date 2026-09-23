function [deltaBest, CBest] = worstCaseOrderedDelta(tau, epsVec, sMin, lambda, M, nGrid)
%WORSTCASEORDEREDDELTA Find delta_i (i in M) in [-eps_i,eps_i] that
%   maximizes the true system cost trueCost(a,tau), where the landing
%   order is taken from the (possibly lied-about) report tauHat = tau,
%   with tauHat(M) = tau(M)+delta, via reorderedSchedule.m.
%
%   This is a genuine "malicious" agent: unlike a fixed delta=-eps (most
%   aggressive under-report), the actual worst delta for the SYSTEM depends
%   on the scenario -- e.g. a vehicle sitting in a tight downstream gap may
%   find that under-reporting *relieves* pressure and *helps* the system,
%   so the true worst case can even be delta>0 for some vehicles.
%
%   Because reorderedSchedule sorts by tauHat, the map delta -> a is
%   piecewise linear (kinks only where the sort order changes) and
%   trueCost(a,tau) is quadratic in a, so the composed objective is
%   piecewise quadratic in delta. A grid search over the box finds the
%   right piece (and is already exact if the order never changes within
%   it, which solveSchedule's QP guarantees is convex there); fmincon then
%   polishes the grid optimum to convergence within that piece.

tau      = tau(:);
epsVec   = epsVec(:);
M        = M(:);
nM       = numel(M);

    function C = negObjective(deltaM)
        tauHat = tau; tauHat(M) = tau(M) + deltaM(:);
        a = reorderedSchedule(tauHat, tau, epsVec, lambda, sMin);
        C = -trueCost(a, tau);
    end

% ---- stage 1: grid search over the box, to locate the right piece ----
grids = cell(nM,1);
for k = 1:nM
    grids{k} = linspace(-epsVec(M(k)), epsVec(M(k)), nGrid);
end
nCombo = nGrid^nM;
comboIdx = zeros(nCombo, nM);
remaining = (0:nCombo-1)';
for k = 1:nM
    comboIdx(:,k) = mod(remaining, nGrid) + 1;
    remaining     = floor(remaining / nGrid);
end

CBestNeg  = inf;
deltaBest = zeros(nM,1);
for c = 1:nCombo
    deltaM = zeros(nM,1);
    for k = 1:nM
        deltaM(k) = grids{k}(comboIdx(c,k));
    end
    Cneg = negObjective(deltaM);
    if Cneg < CBestNeg
        CBestNeg  = Cneg;
        deltaBest = deltaM;
    end
end

% ---- stage 2: polish with fmincon within that piece ----
opts = optimoptions('fmincon', 'Display', 'off');
lb = -epsVec(M);
ub =  epsVec(M);
try
    [deltaPolished, CnegPolished] = fmincon(@negObjective, deltaBest, [], [], [], [], lb, ub, [], opts);
    if CnegPolished <= CBestNeg
        deltaBest = deltaPolished;
        CBestNeg  = CnegPolished;
    end
catch
    % fall back silently on the grid result if fmincon/toolbox unavailable
end

CBest = -CBestNeg;

end
