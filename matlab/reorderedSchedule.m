function [a, order] = reorderedSchedule(tauHat, tauTilde, epsVec, lambda, sMin)
%REORDEREDSCHEDULE Determine the landing order by sorting the reported
%   arrival times (tauHat) ascending, solve the fixed-order robust QP
%   (solveSchedule.m) in that order, then map the result back to the
%   original vehicle indices.
%
%   This is the mechanism by which a vehicle that under-reports its ETA
%   (delta < 0) can try to move earlier in the queue: the order itself is
%   always taken from the raw report (nobody else has any way to order
%   vehicles), for both the robust and non-robust coordinator. What
%   differs is the lambda used inside solveSchedule, which controls how
%   much the resulting assigned time actually reflects that report.

N = numel(tauHat);
tauHat   = tauHat(:);
tauTilde = tauTilde(:);
epsVec   = epsVec(:);
lambda   = lambda(:);

[~, order] = sort(tauHat);

aOrdered = solveSchedule(tauHat(order), tauTilde(order), epsVec(order), lambda(order), sMin);

a = zeros(N,1);
a(order) = aOrdered;

end
