%RUNEPSILONSWEEP Optimality gap (%) AND raw cost vs. lambda-distance,
%   repeated across several fixed surveillance-uncertainty values epsilon.
%
%   Each epsilon value runs the same experiment as runLambdaDistanceSweep.m
%   (sigma-sweep noise on lambdaHat, equalized distance bins over [0,1.2])
%   with scenario.eps held at that single uniform value. All five values
%   are computed fresh this run (raw achieved/oracle cost is now tracked
%   in addition to the gap percentage, so the earlier eps=0.75 result
%   cannot be reused as-is).
%
%   Produces two figures:
%     1) gap (%) vs distance, one line per epsilon
%     2) raw cost vs distance, one pair of lines per epsilon
%        (solid = achieved C(lambdaHat), dashed = oracle C*, same color)

clear; clc;
here = fileparts(mfilename('fullpath'));

epsList = [0.25 0.5 0.75 1.0 1.25];
sigmaList = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
nTrialsPerSigma = 6000;

results = cell(numel(epsList),1);
for i = 1:numel(epsList)
    epsVal = epsList(i);
    fprintf('eps=%.2f: running sweep\n', epsVal);
    seed = round(epsVal*1000) + 100;
    results{i} = runDistanceSweepCore(epsVal, sigmaList, nTrialsPerSigma, seed);
end

colors = lines(numel(epsList));

% ---- Figure 1: gap (%) vs distance ----
fig1 = figure; hold on; grid on;
for i = 1:numel(epsList)
    r = results{i};
    v = ~isnan(r.meanGapPct);
    errorbar(r.binCenter(v), r.meanGapPct(v), r.semGapPct(v), ...
        '-o', 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('\\epsilon = %.2f (n=%d)', epsList(i), r.nEqual));
end
xlabel('||\lambda hat - \lambda^{true}||_2');
ylabel('Optimality gap  (%)');
title('Optimality gap (%) vs. \lambda distance, across surveillance uncertainty \epsilon  (error bars = \pm1 SEM)');
legend('Location', 'northwest');
exportgraphics(fig1, fullfile(here, 'epsilonSweepResults.png'));
exportgraphics(fig1, fullfile(here, 'epsilonSweepResults.pdf'), 'ContentType', 'vector');

% ---- Figure 2: raw cost vs distance (achieved vs oracle) ----
fig2 = figure; hold on; grid on;
for i = 1:numel(epsList)
    r = results{i};
    v = ~isnan(r.meanChat);
    errorbar(r.binCenter(v), r.meanChat(v), r.semChat(v), ...
        '-o', 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('C(\\lambda hat), \\epsilon=%.2f', epsList(i)));
    errorbar(r.binCenter(v), r.meanCoracle(v), r.semCoracle(v), ...
        '--s', 'LineWidth', 1.2, 'Color', colors(i,:), ...
        'HandleVisibility', 'off');
end
xlabel('||\lambda hat - \lambda^{true}||_2');
ylabel('Raw true cost   \Sigma_i (a_i-\tau_i)^2');
title('Raw achieved cost C(\lambda hat) (solid) vs. oracle C^* (dashed, unlabeled), by \epsilon  (error bars = \pm1 SEM)');
legend('Location', 'northwest');
exportgraphics(fig2, fullfile(here, 'epsilonSweepRawCost.png'));
exportgraphics(fig2, fullfile(here, 'epsilonSweepRawCost.pdf'), 'ContentType', 'vector');

save(fullfile(here, 'epsilonSweepResults.mat'), 'epsList', 'results', 'sigmaList', 'nTrialsPerSigma');
