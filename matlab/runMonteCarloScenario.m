%RUNMONTECARLOSCENARIO Monte Carlo over randomized scenarios: instead of
%   the fixed prototype tau=[10 12 13 17 19], eps=0.75, each trial draws a
%   fresh scenario and compares the same 4 coordinators x 3 behaviors as
%   runImpreciseLambdaComparison.m:
%
%     tau_i     ~ Uniform(tauLow, tauHigh), iid per vehicle, NOT sorted
%                 (vehicle identity is just its draw index; M=[2,4] means
%                 "the 2nd and 4th vehicle drawn", whatever tau they land on)
%     eps       ~ Uniform(epsLow, epsHigh), ONE shared draw per trial
%                 (same eps for all 5 vehicles that trial, fresh every trial)
%     sMin      fixed at 2, as in the prototype scenario
%
%   Coordinators (same lambda vectors as runImpreciseLambdaComparison.m):
%     oracle robust, imprecise robust (err 0.3), imprecise robust (err 0.6,
%     crosses the 0.5 dead-zone threshold), non-robust.
%   Behaviors: malicious (worst-case delta_M via worstCaseOrderedDelta.m,
%   found fresh per trial per coordinator since it depends on tau/eps),
%   random (one Uniform(-eps,eps) draw per trial), honest (delta=0).
%
%   Every trial also records ratio = var(tau)/eps -- "how spread out the
%   arrivals are relative to the surveillance resolution" -- for the
%   scatter-vs-ratio figure.
%
%   Output:
%     Figure 1: grouped boxchart + swarmchart of system cost, x=behavior,
%               grouped/colored by coordinator (all trials pooled)
%     Figure 2: system cost vs. ratio=var(tau)/eps, one scatter subplot per
%               behavior, colored by coordinator

clear; clc;
N    = 5;
sMin = 2;
M    = [2; 4];

lambdaOracle     = [1; 0; 1; 0; 1];
lambdaImprecise1 = [0.7; 0.3; 0.7; 0.3; 0.7];
lambdaImprecise2 = [0.4; 0.6; 0.4; 0.6; 0.4];
lambdaNonRobust  = ones(N,1);

baseNames = {'oracle robust', 'imprecise robust (err 0.3)', 'imprecise robust (err 0.6)', 'non-robust'};
lambdas   = {lambdaOracle, lambdaImprecise1, lambdaImprecise2, lambdaNonRobust};
nCoord    = numel(lambdas);
coordNames = cell(1,nCoord);
for ci = 1:nCoord
    coordNames{ci} = sprintf('%s (||\\Delta\\lambda||=%.2f)', baseNames{ci}, norm(lambdas{ci}-lambdaOracle));
end

behaviorNames = {'malicious', 'random', 'honest'};

tauLow = 0; tauHigh = 10;   % expected adjacent-gap ~ tauHigh/6 ~ 1.67 < sMin=2,
                            % so separation binds in most trials while still
                            % leaving some looser scenarios for ratio variety
epsLow = 0.25; epsHigh = 1.25;
nGrid  = 15;
nTrials = 1000;

rng(42);

sysCost = nan(nTrials, nCoord, 3);
ratio   = nan(nTrials, 1);

t0 = tic;
for trial = 1:nTrials
    tau    = tauLow + (tauHigh-tauLow)*rand(N,1);
    epsVal = epsLow + (epsHigh-epsLow)*rand();
    epsVec = epsVal * ones(N,1);
    ratio(trial) = var(tau) / epsVal;

    randDelta = (2*rand(numel(M),1)-1) .* epsVec(M);

    for ci = 1:nCoord
        lambda = lambdas{ci};
        deltaWorst = worstCaseOrderedDelta(tau, epsVec, sMin, lambda, M, nGrid);

        for bi = 1:3
            switch behaviorNames{bi}
                case 'malicious', deltaM = deltaWorst;
                case 'random',    deltaM = randDelta;
                case 'honest',    deltaM = zeros(numel(M),1);
            end
            tauHat = tau; tauHat(M) = tau(M) + deltaM;
            a = reorderedSchedule(tauHat, tau, epsVec, lambda, sMin);
            sysCost(trial, ci, bi) = trueCost(a, tau);
        end
    end

    if mod(trial, 100) == 0
        fprintf('trial %d/%d done (%.1f sec elapsed)\n', trial, nTrials, toc(t0));
    end
end
fprintf('Total time: %.1f sec\n', toc(t0));

here = fileparts(mfilename('fullpath'));
save(fullfile(here, 'monteCarloScenarioResults.mat'), ...
    'sysCost', 'ratio', 'coordNames', 'behaviorNames', 'lambdas', 'M', 'sMin', ...
    'tauLow', 'tauHigh', 'epsLow', 'epsHigh', 'nTrials', 'nGrid');

% ---- Figure 1: boxplot + scatter, x=behavior, group=coordinator ----
%      Boxes and their scatter points are placed at explicit, matching
%      numeric x positions per (behavior, coordinator) pair -- swarmchart's
%      'GroupByColor' isn't supported in this MATLAB version, so grouping
%      is done manually instead of relying on boxchart+swarmchart
%      auto-grouping.
xCenters     = 1:3;
groupOffsets = linspace(-0.3, 0.3, nCoord);
boxWidth     = 0.16;
jitterHalf   = 0.06;
colors = lines(nCoord);

fig1 = figure; hold on; grid on;
boxHandles = gobjects(nCoord,1);
for ci = 1:nCoord
    for bi = 1:3
        x = xCenters(bi) + groupOffsets(ci);
        v = sysCost(:,ci,bi);
        h = boxchart(repmat(x, nTrials, 1), v, ...
            'BoxFaceColor', colors(ci,:), 'BoxWidth', boxWidth, 'MarkerStyle', 'none');
        if bi == 1
            boxHandles(ci) = h;
        end
        jitterX = x + (rand(nTrials,1)-0.5)*2*jitterHalf;
        scatter(jitterX, v, 6, colors(ci,:), 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    end
end
set(gca, 'XTick', xCenters, 'XTickLabel', behaviorNames);
xlim([1-0.5, 3+0.5]);
legend(boxHandles, coordNames, 'Location', 'northeast', 'Interpreter', 'tex');
ylabel('True system cost  \Sigma_i (a_i-\tau_i)^2');
title(sprintf('Monte Carlo (%d random scenarios): system cost by behavior and coordinator', nTrials));

exportgraphics(fig1, fullfile(here, 'monteCarloBoxplot.png'));
exportgraphics(fig1, fullfile(here, 'monteCarloBoxplot.pdf'), 'ContentType', 'vector');

% ---- Figure 2: system cost vs. ratio=var(tau)/eps, one subplot per
%      behavior, coordinators overlaid by color. Y-axis unified across
%      subplots; legend uses short coordinator names (no ||Delta lambda||,
%      that's already in the companion boxplot/table) ----
shortNames = regexprep(coordNames, '\s*\(\|\|\\Delta\\lambda\|\|=[0-9.]+\)', '');
fig2 = figure;
colors = lines(nCoord);
yMax = max(sysCost(:), [], 'omitnan') * 1.05;
axHandles = gobjects(1,3);
for bi = 1:3
    axHandles(bi) = subplot(1,3,bi); hold on; grid on;
    for ci = 1:nCoord
        scatter(ratio, sysCost(:,ci,bi), 8, colors(ci,:), 'filled', 'MarkerFaceAlpha', 0.35);
    end
    ylim([0 yMax]);
    xlabel('var(\tau) / \epsilon');
    ylabel('True system cost');
    title(behaviorNames{bi});
end
linkaxes(axHandles, 'xy');
subplot(1,3,1);
legend(shortNames, 'Location', 'northeast');
sgtitle(sprintf('System cost vs. arrival spread / surveillance resolution (%d trials)', nTrials));

exportgraphics(fig2, fullfile(here, 'monteCarloScatterVsRatio.png'));
exportgraphics(fig2, fullfile(here, 'monteCarloScatterVsRatio.pdf'), 'ContentType', 'vector');
