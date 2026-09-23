%REPLOTMONTECARLOSCATTER Regenerate the ratio-scatter figure (overlaid,
%   one subplot per behavior, all coordinators overlaid by color) from the
%   already-saved runMonteCarloScenario.m results. Y-axis is unified
%   across the 3 subplots, and the legend uses short coordinator names
%   (the ||Delta lambda|| distance is dropped -- it's already in the
%   companion boxplot/table). Keep this in sync with the Figure 2 block in
%   runMonteCarloScenario.m.

clear; clc;
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'monteCarloScenarioResults.mat'));
sysCost = S.sysCost; ratio = S.ratio; coordNames = S.coordNames;
behaviorNames = S.behaviorNames; nTrials = S.nTrials;
nCoord = numel(coordNames);

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
