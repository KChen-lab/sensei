clear

a = [3 4 5 6]';

p = gamrnd(repmat(a, [1 10000000]), 1);
p = p ./ sum(p);

p0 = p(2, :) + p(4, :);
p1 = p(3, :) + p(4, :);


empirical_corr = corr([p0; p1]')
theoretical_corr = (a(1) * a(4) - a(2) * a(3)) / ...
    sqrt((a(1) + a(2)) * (a(1) + a(3)) * (a(4) + a(2)) * (a(4) + a(3)))