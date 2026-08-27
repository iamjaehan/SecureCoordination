function c = trueCost(a, tau)
%TRUECOST System-level cost evaluated against ground-truth ETAs, sum_i (a_i-tau_i)^2.
c = sum((a(:) - tau(:)).^2);
end
