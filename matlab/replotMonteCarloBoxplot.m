%REPLOTMONTECARLOBOXPLOT Regenerate the grouped boxplot+scatter figure from
%   the already-saved runMonteCarloScenario.m results, without rerunning
%   the (slow) simulation. Keep this in sync with the Figure 1 block in
%   runMonteCarloScenario.m.
%
%   Boxes and their scatter points are placed at explicit, matching
%   numeric x positions per (behavior, coordinator) pair -- swarmchart's
%   'GroupByColor' isn't supported in this MATLAB version, so this is done
%   manually instead of relying on boxchart+swarmchart auto-grouping.

clear; clc;
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'monteCarloScenarioResults.mat'));
sysCost = S.sysCost; coordNames = S.coordNames; behaviorNames = S.behaviorNames;
nTrials = S.nTrials;
nCoord = numel(coordNames);

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
