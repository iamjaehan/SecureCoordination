%REPLOTMONTECARLONORMALIZED Same grouped box+scatter layout as
%   replotMonteCarloBoxplot.m, but each trial's cost is normalized by that
%   SAME trial's oracle-robust cost under the SAME behavior:
%
%     normCost(trial,ci,bi) = sysCost(trial,ci,bi) / sysCost(trial,oracle,bi)
%
%   This cancels out trial-to-trial difficulty (some random scenarios are
%   just inherently harder/easier than others) and isolates each
%   coordinator's relative performance against the oracle's own defense
%   level for that exact trial. Oracle itself is therefore always exactly
%   1. Trials where the oracle's cost is ~0 for a given behavior (no
%   separation constraint bound at all) are dropped for that behavior,
%   since the ratio is undefined/meaningless there.

clear; clc;
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'monteCarloScenarioResults.mat'));
sysCost = S.sysCost; coordNames = S.coordNames; behaviorNames = S.behaviorNames;
nTrials = S.nTrials;
nCoord = numel(coordNames);
oracleIdx = 1;   % lambdas = {lambdaOracle, lambdaImprecise1, lambdaImprecise2, lambdaNonRobust}

zeroTol = 1e-6;
normCost = nan(nTrials, nCoord, 3);
nDropped = zeros(1,3);
for bi = 1:3
    oracleCost = sysCost(:, oracleIdx, bi);
    valid = oracleCost > zeroTol;
    nDropped(bi) = nTrials - nnz(valid);
    for ci = 1:nCoord
        normCost(valid, ci, bi) = sysCost(valid, ci, bi) ./ oracleCost(valid);
    end
end
fprintf('Trials dropped per behavior (oracle cost ~0): malicious=%d, random=%d, honest=%d (of %d)\n', ...
    nDropped(1), nDropped(2), nDropped(3), nTrials);

xCenters     = 1:3;
groupOffsets = linspace(-0.3, 0.3, nCoord);
boxWidth     = 0.10;
jitterHalf   = 0.035;
colors = lines(nCoord);

fig3 = figure; hold on; grid on;
set(gca, 'YScale', 'log');
yline(1, '--k', 'oracle = 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
boxHandles = gobjects(nCoord,1);
for ci = 1:nCoord
    for bi = 1:3
        x = xCenters(bi) + groupOffsets(ci);
        v = normCost(:,ci,bi);
        v = v(~isnan(v));
        h = boxchart(repmat(x, numel(v), 1), v, ...
            'BoxFaceColor', colors(ci,:), 'BoxWidth', boxWidth, 'MarkerStyle', 'none');
        if bi == 1
            boxHandles(ci) = h;
        end
        jitterX = x + (rand(numel(v),1)-0.5)*2*jitterHalf;
        scatter(jitterX, v, 6, colors(ci,:), 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    end
end
set(gca, 'XTick', xCenters, 'XTickLabel', behaviorNames);
xlim([1-0.5, 3+0.5]);
legend(boxHandles, coordNames, 'Location', 'northeast', 'Interpreter', 'tex');
ylabel('system cost / oracle-robust cost  (same trial & behavior)');
title(sprintf('Monte Carlo (%d trials): cost normalized by oracle-robust, per trial per behavior', nTrials));
ylim([0.6 10]);

exportgraphics(fig3, fullfile(here, 'monteCarloNormalizedBoxplot.png'));
exportgraphics(fig3, fullfile(here, 'monteCarloNormalizedBoxplot.pdf'), 'ContentType', 'vector');
