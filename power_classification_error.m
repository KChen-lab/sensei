M = [10, 10];
N = [1000, 1000];

a = [7, 10];
b = [10, 7];

alpha = 0.05;

q_FN = 0.2;
q_FP = 0.1;
%figure;
%hold on;
%x = 0:.01:1;
%for i = 1:2
%    plot(x, betapdf(x, a(i), b(i)))
%end

Ep = (1 - q_FP - q_FN) * a ./ (a + b) + q_FP;
Vp = (1 - q_FP - q_FN) ^ 2 * a .* b .* (a + b + N) ./ ((N .* (a + b) .^ 2) .* (a + b + 1));

S = sqrt(sum((M - 1) .* Vp) ./ (sum(M) - 2));

Et = diff(Ep) / S / sqrt(sum(1 ./ M));

t_star = tinv(1 - alpha ./ [1 2], sum(M) - 2);

beta = normcdf(t_star - Et); % false negative rate for 1-sided and 2-sided t-tests

beta