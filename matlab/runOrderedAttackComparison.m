%RUNORDEREDATTACKCOMPARISON Robust vs. non-robust coordination, compared
%   across three behavior models for a fixed subset M of vehicles capable
%   of misreporting their ETA to move earlier in the landing queue:
%
%     - "malicious": delta_M in [-eps,eps]^|M| chosen to MAXIMIZE the true
%                    system cost (via worstCaseOrderedDelta.m: grid search
%                    + fmincon polish), found separately per coordinator
%                    since the worst delta depends on lambda
%     - "random"   : delta_i ~ Uniform(-eps_i, eps_i) for i in M
%     - "honest"   : delta_i = 0 for i in M
%
%   In both coordinators, the landing ORDER is taken from the raw report
%   tauHat (ascending sort) -- there is no other way to order vehicles.
%   What differs is the lambda used by solveSchedule once that order is
%   fixed:
%     - robust:     lambda(i)=0 for i in M (oracle distrust), 1 otherwise
%     - non-robust: lambda(i)=1 for everyone (fully trust every report)
%
%   Since robust ignores tauHat_i entirely for i in M, its resulting cost
%   for those vehicles is driven purely by tauTilde_i and eps_i -- constant
%   regardless of M's actual behavior, as long as the sort order doesn't
%   change. Non-robust instead tracks tauHat_i directly, so its cost grows
%   with how aggressively M lies. This script quantifies where (if
%   anywhere) that crossover happens for this scenario.

clear; clc;
scenario = defaultScenario();
N        = scenario.N;
tau      = scenario.tau;
epsVec   = scenario.eps;
sMin     = scenario.sMin;

M = [2; 4];   % vehicles capable of misreporting; coordinator knows this set

lambdaOracle    = ones(N,1); lambdaOracle(M) = 0;   % robust: distrust M
lambdaNonRobust = ones(N,1);                         % non-robust: trust everyone
tauTilde        = tau;                               % exact surveillance here

coordNames = {'robust', 'non-robust'};
lambdas    = {lambdaOracle, lambdaNonRobust};
behaviorNames = {'malicious', 'random', 'honest'};

nTrials = 5000;
nGrid   = 41;
rng(11);

sysCostMean = nan(2,3); sysCostSem = zeros(2,3);
ownCostMean = nan(2,3); ownCostSem = zeros(2,3);
orderEverChanged = false(2,3);
worstDelta = cell(2,1);

for ci = 1:2
    lambda = lambdas{ci};
    worstDelta{ci} = worstCaseOrderedDelta(tau, epsVec, sMin, lambda, M, nGrid);
    for bi = 1:3
        switch behaviorNames{bi}
            case 'malicious'
                deltaDraws = { worstDelta{ci} };
            case 'honest'
                deltaDraws = { zeros(numel(M),1) };
            case 'random'
                deltaDraws = cell(nTrials,1);
                for t = 1:nTrials
                    deltaDraws{t} = (2*rand(numel(M),1)-1) .* epsVec(M);
                end
        end

        nS = numel(deltaDraws);
        Csys = zeros(nS,1);
        Cown = zeros(nS,1);
        for s = 1:nS
            deltaM = deltaDraws{s};
            tauHat = tau; tauHat(M) = tau(M) + deltaM;

            [a, order] = reorderedSchedule(tauHat, tauTilde, epsVec, lambda, sMin);

            Csys(s) = trueCost(a, tau);
            Cown(s) = mean((a(M) - tau(M)).^2);
            if ~isequal(order, (1:N)')
                orderEverChanged(ci,bi) = true;
            end
        end

        sysCostMean(ci,bi) = mean(Csys); sysCostSem(ci,bi) = std(Csys)/sqrt(nS);
        ownCostMean(ci,bi) = mean(Cown); ownCostSem(ci,bi) = std(Cown)/sqrt(nS);
    end
end

fprintf('%-12s %-12s %10s %12s %10s\n', 'coordinator', 'behavior', 'sysCost', 'ownCost(M)', 'orderChg?');
for ci = 1:2
    for bi = 1:3
        fprintf('%-12s %-12s %10.4f %12.4f %10s\n', coordNames{ci}, behaviorNames{bi}, ...
            sysCostMean(ci,bi), ownCostMean(ci,bi), mat2str(orderEverChanged(ci,bi)));
    end
    fprintf('  -> worst-case delta_M for %s: [%s]\n', coordNames{ci}, num2str(worstDelta{ci}'));
end

% ---- Figure 1: system cost, grouped by behavior ----
fig1 = figure; hold on; grid on;
b = bar(sysCostMean');
nbars = numel(b);
xpos = nan(3,nbars);
for k = 1:nbars
    xpos(:,k) = b(k).XEndPoints;
end
errorbar(xpos, sysCostMean', sysCostSem', 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:3, 'XTickLabel', behaviorNames);
legend(coordNames, 'Location', 'northwest');
ylabel('True system cost  \Sigma_i (a_i-\tau_i)^2');
title(sprintf('System cost: robust vs. non-robust, M=[%s]', num2str(M')));

here = fileparts(mfilename('fullpath'));
exportgraphics(fig1, fullfile(here, 'orderedAttackSystemCost.png'));
exportgraphics(fig1, fullfile(here, 'orderedAttackSystemCost.pdf'), 'ContentType', 'vector');

% ---- Figure 2: the misreporting vehicles' own cost ----
fig2 = figure; hold on; grid on;
b2 = bar(ownCostMean');
xpos2 = nan(3,nbars);
for k = 1:nbars
    xpos2(:,k) = b2(k).XEndPoints;
end
errorbar(xpos2, ownCostMean', ownCostSem', 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:3, 'XTickLabel', behaviorNames);
legend(coordNames, 'Location', 'northwest');
ylabel('Mean own cost of M,  mean_i\in M (a_i-\tau_i)^2');
title(sprintf('Misreporting vehicles'' own outcome, M=[%s]', num2str(M')));

exportgraphics(fig2, fullfile(here, 'orderedAttackOwnCost.png'));
exportgraphics(fig2, fullfile(here, 'orderedAttackOwnCost.pdf'), 'ContentType', 'vector');

save(fullfile(here, 'orderedAttackComparisonResults.mat'), ...
    'M', 'lambdaOracle', 'lambdaNonRobust', 'nTrials', 'worstDelta', ...
    'sysCostMean', 'sysCostSem', 'ownCostMean', 'ownCostSem', 'orderEverChanged', ...
    'coordNames', 'behaviorNames');
