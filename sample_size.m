M = [5, 5];
N = [1000, 800];

a = [8, 10];
b = [10, 8];

alpha = 0.05;

%figure;
%hold on;
%x = 0:.01:1;
%for i = 1:2
%    plot(x, betapdf(x, a(i), b(i)))
%end

Ep = a ./ (a + b);
Vp = a .* b .* (a + b + N) ./ ((M .* N .* (a + b) .^ 2) .* (a + b + 1));

S = sqrt(sum((M - 1) .* Vp) ./ (sum(M) - 2));

Et = diff(Ep) / S / sqrt(sum(1 ./ M));

t_star = tinv(1 - alpha ./ [1 2], sum(M) - 2);

beta = normcdf(t_star - Et); % false negative rate for 1-sided and 2-sided t-tests

beta