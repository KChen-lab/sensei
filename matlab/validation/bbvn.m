clear;
close all;

% Parameters
a = 1;
b = 4;
N = 100;

x = 0:N;

%% Normal distr
Ep = a / (a + b);
Vp = a * b * (a + b + N) / N / (a + b) ^ 2 / (a + b + 1);

n_cdf = normcdf(x / N, Ep, sqrt(Vp));

%% BetaBinom
bb_pmf = exp(gammaln(N + 1)-gammaln(x + 1)-gammaln(N - x + 1)) .* ...
    beta((a + x),(b + N - x)) / beta(a, b);
bb_cdf = cumsum(bb_pmf);

figure;
plot(x / N, n_cdf)
hold on
stairs(x / N, bb_cdf)

%% Distance measurement
figure, hold on;
N_list = {10 50 100 500 1000};
for ii = 1:length(N_list)
    N = N_list{ii};
    Ep = a / (a + b);
    Vp = a * b * (a + b + N) / N / (a + b) ^ 2 / (a + b + 1);
    
    step = 10000;
    x2 = linspace(-0.5 * N, 1.5 * N, step);
    n_cdf = normcdf(x2 / N, Ep, sqrt(Vp));

    floor_x2 = floor(x2);
    floor_x2(x2 < 0) = 0;
    floor_x2(x2 > N) = N;
    bb_pmf = exp(gammaln(N + 1)-gammaln(floor_x2 + 1)-gammaln(N - floor_x2 + 1)) .* ...
        beta((a + floor_x2),(b + N - floor_x2)) / beta(a, b);
    bb_pmf(x2 < 0) = 0;
    bb_pmf([false bb_pmf(1:end-1) == bb_pmf(2:end)]) = 0;
    bb_cdf = cumsum(bb_pmf);

    dist = sum(abs(n_cdf - bb_cdf)) / (step - 1);

    ax.ColorOrderIndex = ii;
    plot(x2 / N, n_cdf, '--')
    hold on
    ax.ColorOrderIndex = ii;
    stairs(x2 / N, bb_cdf, "-")
    ylim([0 1])
    xlim([-0.5 1.5])
end
legend("10", "50", "100", "500", "1000", "5000");