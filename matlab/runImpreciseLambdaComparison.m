%RUNIMPRECISELAMBDACOMPARISON How much does robust coordination's benefit
%   degrade as the coordinator's trust estimate lambda drifts away from the
%   oracle (i.e. is estimated, not the true 0/1 identity of who is
%   malicious)? Each coordinator is labeled by ||lambda - lambdaOracle||_2
%   so the drift is directly comparable across coordinators, non-robust
%   included (non-robust is just the most extreme point on this axis).
%
%   Coordinators, same landing-order mechanism as runOrderedAttackComparison.m
%   (order = ascending sort of tauHat, then solveSchedule in that order):
%     - oracle robust:          lambda = [1 0 1 0 1]   (perfectly knows M=[2,4])
%     - imprecise robust:       lambda = [0.7 0.3 0.7 0.3 0.7]  (error 0.3,
%                                still below the 0.5 dead-zone threshold --
%                                see the lambda-deadzone finding: worst-case
%                                cost should still match the oracle exactly)
%     - imprecise robust (large err): lambda = [0.4 0.6 0.4 0.6 0.4]
%                                (error 0.6, *crosses* the 0.5 threshold, so
%                                M is now trusted MORE than the honest
%                                vehicles -- worst-case cost should start
%                                to degrade past the oracle's ceiling here)
%     - non-robust:             lambda = [1 1 1 1 1]   (fully trusts every report)
%
%   Same three behavior models for M = [2,4] as runOrderedAttackComparison.m,
%   always compared side by side:
%     - malicious: delta_M chosen to MAXIMIZE true system cost (grid search
%                  + fmincon polish via worstCaseOrderedDelta.m), found
%                  separately for each lambda vector
%     - random:    delta_i ~ Uniform(-eps_i, eps_i), averaged over nTrials
%     - honest:    delta_i = 0 (fixed)

clear; clc;
scenario = defaultScenario();
N        = scenario.N;
tau      = scenario.tau;
epsVec   = scenario.eps;
sMin     = scenario.sMin;

M = [2; 4];

lambdaOracle      = ones(N,1); lambdaOracle(M) = 0;
lambdaImprecise1  = [0.7; 0.3; 0.7; 0.3; 0.7];
lambdaImprecise2  = [0.4; 0.6; 0.4; 0.6; 0.4];
lambdaNonRobust   = ones(N,1);
tauTilde          = tau;

baseNames = {'oracle robust', 'imprecise robust', 'imprecise robust (large err)', 'non-robust'};
lambdas   = {lambdaOracle, lambdaImprecise1, lambdaImprecise2, lambdaNonRobust};
nCoord    = numel(lambdas);

coordNames = cell(1,nCoord);
for ci = 1:nCoord
    lambdaDist = norm(lambdas{ci} - lambdaOracle);
    coordNames{ci} = sprintf('%s (||\\Delta\\lambda||=%.2f)', baseNames{ci}, lambdaDist);
end

behaviorNames = {'malicious', 'random', 'honest'};

nTrials = 5000;
nGrid   = 41;
rng(11);

sysCostMean = nan(nCoord,3); sysCostSem = zeros(nCoord,3);
ownCostMean = nan(nCoord,3); ownCostSem = zeros(nCoord,3);
orderEverChanged = false(nCoord,3);
worstDelta = cell(nCoord,1);

for ci = 1:nCoord
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

fprintf('%-45s %-12s %10s %12s %10s\n', 'coordinator', 'behavior', 'sysCost', 'ownCost(M)', 'orderChg?');
for ci = 1:nCoord
    for bi = 1:3
        fprintf('%-45s %-12s %10.4f %12.4f %10s\n', coordNames{ci}, behaviorNames{bi}, ...
            sysCostMean(ci,bi), ownCostMean(ci,bi), mat2str(orderEverChanged(ci,bi)));
    end
    fprintf('  -> worst-case delta_M: [%s]\n', num2str(worstDelta{ci}'));
end

% ---- Figure 1: system cost ----
fig1 = figure; hold on; grid on;
b = bar(sysCostMean');
nbars = numel(b);
xpos = nan(3,nbars);
for k = 1:nbars
    xpos(:,k) = b(k).XEndPoints;
end
errorbar(xpos, sysCostMean', sysCostSem', 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:3, 'XTickLabel', behaviorNames);
legend(coordNames, 'Location', 'northwest', 'Interpreter', 'tex');
ylabel('True system cost  \Sigma_i (a_i-\tau_i)^2');
title(sprintf('System cost vs. \\lambda precision, M=[%s]', num2str(M')));

here = fileparts(mfilename('fullpath'));
exportgraphics(fig1, fullfile(here, 'impreciseLambdaSystemCost.png'));
exportgraphics(fig1, fullfile(here, 'impreciseLambdaSystemCost.pdf'), 'ContentType', 'vector');

% ---- Figure 2: M's own outcome ----
fig2 = figure; hold on; grid on;
b2 = bar(ownCostMean');
xpos2 = nan(3,nbars);
for k = 1:nbars
    xpos2(:,k) = b2(k).XEndPoints;
end
errorbar(xpos2, ownCostMean', ownCostSem', 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:3, 'XTickLabel', behaviorNames);
legend(coordNames, 'Location', 'northwest', 'Interpreter', 'tex');
ylabel('Mean own cost of M,  mean_i\in M (a_i-\tau_i)^2');
title(sprintf('Misreporting vehicles'' own outcome vs. \\lambda precision, M=[%s]', num2str(M')));

exportgraphics(fig2, fullfile(here, 'impreciseLambdaOwnCost.png'));
exportgraphics(fig2, fullfile(here, 'impreciseLambdaOwnCost.pdf'), 'ContentType', 'vector');

save(fullfile(here, 'impreciseLambdaComparisonResults.mat'), ...
    'M', 'lambdaOracle', 'lambdaImprecise1', 'lambdaImprecise2', 'lambdaNonRobust', 'nTrials', 'worstDelta', ...
    'sysCostMean', 'sysCostSem', 'ownCostMean', 'ownCostSem', 'orderEverChanged', ...
    'coordNames', 'behaviorNames');
