%RUNLAMBDADISTANCESWEEP Optimality gap (%) vs. ||lambdaHat - lambdaTrue||_2.
%
%   C* is the oracle cost: solve Eq. (12) using the *realized* per-round
%   truthful/liar label for each vehicle (lambda_i = 1 if it told the
%   truth this round, 0 if it lied). This oracle is not achievable by any
%   real coordinator, since lying is stochastic (governed by
%   scenario.lambdaTrue) and only the *probability* of lying is ever
%   knowable in advance -- never this round's realized outcome.
%
%   The coordinator's estimate is
%       lambdaHat = clip(lambdaTrue + sigma*randn(N,1), 0, 1),
%   i.e. i.i.d. noise with THE SAME sigma applied to every vehicle (no
%   vehicle is treated differently). sigma itself is swept over sigmaList
%   purely as a mechanism to generate a spread of realized 2-norm
%   distances ||lambdaHat-lambdaTrue||_2 -- sweeping sigma per se is a
%   separate question, left for later.
%
%   Gap is reported as a percentage of the oracle cost, averaged per
%   trial: mean_trial[ 100*(Chat-Coracle)/Coracle ].
%
%   Sample-size control: all (sigma, trial) draws are pooled and binned by
%   realized distance into 0.1-wide bins over [0, 1.2]. Each bin is then
%   randomly DOWNSAMPLED to the same common count nEqual, so every point
%   on the resulting curve is computed from an equal number of samples.

clear; clc;

scenario = defaultScenario();
sigmaList = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
nTrialsPerSigma = 6000;
rng(5);

N = scenario.N;
total = numel(sigmaList) * nTrialsPerSigma;
allDist   = zeros(total,1);
allGapPct = zeros(total,1);
n = 0;

for sigma = sigmaList
    for trial = 1:nTrialsPerSigma
        [tauHat, tauTilde, liarMask] = drawRealization(scenario);

        lambdaOracle = double(~liarMask);
        aOracle = solveSchedule(tauHat, tauTilde, scenario.eps, lambdaOracle, scenario.sMin);
        Coracle = trueCost(aOracle, scenario.tau);

        lambdaHat = min(max(scenario.lambdaTrue + sigma*randn(N,1), 0), 1);
        aHat = solveSchedule(tauHat, tauTilde, scenario.eps, lambdaHat, scenario.sMin);
        Chat = trueCost(aHat, scenario.tau);

        n = n+1;
        allDist(n)   = norm(lambdaHat - scenario.lambdaTrue, 2);
        allGapPct(n) = 100 * (Chat - Coracle) / Coracle;
    end
    fprintf('sigma = %.2f done (n = %d so far)\n', sigma, n);
end

% bin by realized 2-norm distance
edgeMax = 1.2;
edges = 0:0.1:edgeMax;
nBins = numel(edges)-1;
binIdx = discretize(allDist, edges);
binCenter = ((edges(1:end-1)+edges(2:end))/2)';

rawCount = zeros(nBins,1);
for b = 1:nBins
    rawCount(b) = sum(binIdx==b);
end
fprintf('\nraw bin counts (before equalizing):\n');
for b = 1:nBins
    fprintf('  [%.1f,%.1f): n=%d\n', edges(b), edges(b+1), rawCount(b));
end

nEqual = min(rawCount(rawCount > 0));
fprintf('\nequalizing every bin to nEqual = %d samples\n\n', nEqual);

meanGapPct = nan(nBins,1);
stdGapPct  = nan(nBins,1);
semGapPct  = nan(nBins,1);
rng(6); % separate stream for the downsampling draw, for reproducibility
for b = 1:nBins
    idx = find(binIdx==b);
    if numel(idx) >= nEqual && nEqual > 0
        sel = idx(randperm(numel(idx), nEqual));
        meanGapPct(b) = mean(allGapPct(sel));
        stdGapPct(b)  = std(allGapPct(sel));
        semGapPct(b)  = stdGapPct(b) / sqrt(nEqual);
    end
end
valid = ~isnan(meanGapPct);

fig = figure; hold on; grid on;
errorbar(binCenter(valid), meanGapPct(valid), semGapPct(valid), '-o', 'LineWidth', 1.5);
xlabel('||\lambda hat - \lambda^{true}||_2');
ylabel('Optimality gap  (%)');
title(sprintf('Optimality gap (%%) vs. distance from \\lambda^{true}  (n=%d per bin, error bars = \\pm1 SEM)', nEqual));

here = fileparts(mfilename('fullpath'));
exportgraphics(fig, fullfile(here, 'lambdaDistanceSweepResults.png'));
exportgraphics(fig, fullfile(here, 'lambdaDistanceSweepResults.pdf'), 'ContentType', 'vector');
save(fullfile(here, 'lambdaDistanceSweepResults.mat'), ...
    'sigmaList', 'nTrialsPerSigma', 'allDist', 'allGapPct', 'edges', ...
    'binCenter', 'rawCount', 'nEqual', 'meanGapPct', 'stdGapPct', 'semGapPct', 'scenario');
