function res = runDistanceSweepCore(epsVal, sigmaList, nTrialsPerSigma, seed)
%RUNDISTANCESWEEPCORE Optimality gap (%) AND raw cost vs. ||lambdaHat-lambdaTrue||_2,
%   for a fixed surveillance uncertainty epsVal applied uniformly to all
%   vehicles. sigma is swept only to generate a spread of realized
%   distances (same sigma applied to every vehicle each draw); every
%   0.1-wide distance bin over [0,1.2] is downsampled to a common count
%   nEqual so all bins carry equal sample size.
%
%   res fields:
%     binCenter               distance-bin centers
%     meanGapPct, stdGapPct, semGapPct   100*(Chat-Coracle)/Coracle stats
%     meanChat,  stdChat,  semChat       raw achieved cost C(lambdaHat)
%     meanCoracle, stdCoracle, semCoracle  raw oracle cost C*
%     nEqual, rawCount

scenario = defaultScenario();
scenario.eps = epsVal * ones(scenario.N,1);
N = scenario.N;

rng(seed);
total = numel(sigmaList) * nTrialsPerSigma;
allDist    = zeros(total,1);
allGapPct  = zeros(total,1);
allChat    = zeros(total,1);
allCoracle = zeros(total,1);
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
        allDist(n)    = norm(lambdaHat - scenario.lambdaTrue, 2);
        allGapPct(n)  = 100 * (Chat - Coracle) / Coracle;
        allChat(n)    = Chat;
        allCoracle(n) = Coracle;
    end
    fprintf('  eps=%.2f: sigma=%.2f done (n=%d)\n', epsVal, sigma, n);
end

edges = 0:0.1:1.2;
nBins = numel(edges)-1;
binIdx = discretize(allDist, edges);
binCenter = ((edges(1:end-1)+edges(2:end))/2)';

rawCount = zeros(nBins,1);
for b = 1:nBins
    rawCount(b) = sum(binIdx==b);
end
nEqual = min(rawCount(rawCount > 0));

meanGapPct = nan(nBins,1); stdGapPct = nan(nBins,1); semGapPct = nan(nBins,1);
meanChat   = nan(nBins,1); stdChat   = nan(nBins,1); semChat   = nan(nBins,1);
meanCoracle = nan(nBins,1); stdCoracle = nan(nBins,1); semCoracle = nan(nBins,1);

rng(seed+1);
for b = 1:nBins
    idx = find(binIdx==b);
    if numel(idx) >= nEqual && nEqual > 0
        sel = idx(randperm(numel(idx), nEqual));

        meanGapPct(b) = mean(allGapPct(sel));
        stdGapPct(b)  = std(allGapPct(sel));
        semGapPct(b)  = stdGapPct(b) / sqrt(nEqual);

        meanChat(b) = mean(allChat(sel));
        stdChat(b)  = std(allChat(sel));
        semChat(b)  = stdChat(b) / sqrt(nEqual);

        meanCoracle(b) = mean(allCoracle(sel));
        stdCoracle(b)  = std(allCoracle(sel));
        semCoracle(b)  = stdCoracle(b) / sqrt(nEqual);
    end
end

res.binCenter = binCenter;
res.nEqual = nEqual;
res.rawCount = rawCount;
res.meanGapPct = meanGapPct; res.stdGapPct = stdGapPct; res.semGapPct = semGapPct;
res.meanChat = meanChat; res.stdChat = stdChat; res.semChat = semChat;
res.meanCoracle = meanCoracle; res.stdCoracle = stdCoracle; res.semCoracle = semCoracle;

end
