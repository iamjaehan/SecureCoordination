%RUNLAMBDAERRORSWEEP Optimality gap vs. lambda-estimation error.
%
%   Oracle case: lambda_i = 0 for a vehicle that actually lied this round,
%                lambda_i = 1 for a vehicle that reported truthfully this round.
%   Solving Eq. (12) with the oracle lambda gives a* and the optimal
%   true-cost C* = sum_i (a*_i - tau_i)^2.
%
%   Estimation error: each vehicle's oracle label is independently flipped
%   with probability pErr before being handed to the coordinator. Sweeping
%   pErr from 0 (perfect knowledge) upward traces the optimality gap
%   E[C(lambdaHat) - C*].

clear; clc;

scenario = defaultScenario();
pErrList = 0:0.05:0.5;
nTrials  = 2000;
rng(1);

meanGap = zeros(size(pErrList));
semGap  = zeros(size(pErrList));

for k = 1:numel(pErrList)
    pErr = pErrList(k);
    gaps = zeros(nTrials,1);

    for trial = 1:nTrials
        [tauHat, tauTilde, liarMask] = drawRealization(scenario);

        lambdaOracle = double(~liarMask);
        aOracle = solveSchedule(tauHat, tauTilde, scenario.eps, lambdaOracle, scenario.sMin);
        Coracle = trueCost(aOracle, scenario.tau);

        flip = rand(scenario.N,1) < pErr;
        lambdaHat = lambdaOracle;
        lambdaHat(flip) = 1 - lambdaHat(flip);
        aHat = solveSchedule(tauHat, tauTilde, scenario.eps, lambdaHat, scenario.sMin);
        Chat = trueCost(aHat, scenario.tau);

        gaps(trial) = Chat - Coracle;
    end

    meanGap(k) = mean(gaps);
    semGap(k)  = std(gaps) / sqrt(nTrials);
    fprintf('pErr = %.2f : mean gap = %.4f (sem %.4f)\n', pErr, meanGap(k), semGap(k));
end

fig = figure; hold on; grid on;
errorbar(pErrList, meanGap, semGap, '-o', 'LineWidth', 1.5);

% overlay the irreducible statistical-only baseline (lambda = lambda_true,
% no knowledge of this round's realized outcome), if it has been computed
here = fileparts(mfilename('fullpath'));
statFile = fullfile(here, 'statisticalGapResult.mat');
if isfile(statFile)
    s = load(statFile, 'meanGapStat');
    yline(s.meanGapStat, '--', sprintf('\\lambda=\\lambda^{true} only (%.3f)', s.meanGapStat), ...
        'LabelHorizontalAlignment', 'left');
end

xlabel('\lambda misclassification probability p_{err}');
ylabel('Optimality gap  E[C(\lambda hat) - C^*]');
title('Optimality gap vs. \lambda estimation error');

exportgraphics(fig, fullfile(here, 'lambdaSweepResults.png'));
exportgraphics(fig, fullfile(here, 'lambdaSweepResults.pdf'), 'ContentType', 'vector');
save(fullfile(here, 'lambdaSweepResults.mat'), ...
    'pErrList', 'meanGap', 'semGap', 'scenario', 'nTrials');
