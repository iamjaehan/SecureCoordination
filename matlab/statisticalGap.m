%STATISTICALGAP Irreducible gap from knowing only lambda_true (a-priori lie
%   probability) versus the oracle, which knows this round's realized
%   truthful/liar outcome. This isolates the loss due to stochastic lying
%   itself, separate from any lambda-estimation/misclassification error.

clear; clc;
scenario = defaultScenario();
rng(2);
nTrials = 5000;

gapsStat = zeros(nTrials,1);
for trial = 1:nTrials
    [tauHat, tauTilde, liarMask] = drawRealization(scenario);

    lambdaOracle = double(~liarMask);
    aOracle = solveSchedule(tauHat, tauTilde, scenario.eps, lambdaOracle, scenario.sMin);
    Coracle = trueCost(aOracle, scenario.tau);

    % coordinator knows only the a-priori lie probability lambdaTrue,
    % not this round's realized outcome
    aStat = solveSchedule(tauHat, tauTilde, scenario.eps, scenario.lambdaTrue, scenario.sMin);
    Cstat = trueCost(aStat, scenario.tau);

    gapsStat(trial) = Cstat - Coracle;
end

meanGapStat = mean(gapsStat);
semGapStat  = std(gapsStat) / sqrt(nTrials);
fprintf('statistical-only gap: mean = %.4f (sem %.4f)\n', meanGapStat, semGapStat);

save(fullfile(fileparts(mfilename('fullpath')), 'statisticalGapResult.mat'), ...
    'meanGapStat', 'semGapStat', 'nTrials');
