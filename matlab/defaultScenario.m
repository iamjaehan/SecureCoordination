function scenario = defaultScenario()
%DEFAULTSCENARIO Prototype vertiport sequencing scenario (paper Sec. IV-F).
%   N = 5 vehicles, fixed landing order 1->2->3->4->5.

scenario.N          = 5;
scenario.tau        = [10 12 13 17 19]';   % true ETAs
% Surveillance uncertainty: held at one common value across all vehicles
% for now (mean of the original per-vehicle prototype values
% [0.3 1.2 0.4 1.5 0.3]). Revisit / vary this later.
scenario.eps        = 0.75 * ones(5,1);
scenario.sMin       = 2;                    % minimum separation
scenario.lambdaTrue = [0.9 0.3 0.9 0.2 0.9]'; % P(vehicle i reports truthfully)
scenario.malDelta   = -3 * ones(scenario.N,1); % malicious reporting deviation

end
