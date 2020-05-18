clear;
close all;

%% Distance measurement
mean_list = 0.1:0.1:0.5;
std_list = 0.05:0.01:0.3;
N_list = [10 50 100 500];

res = nan(length(mean_list), length(std_list), length(N_list));

for kk = 1:length(mean_list)
    for jj = 1:length(std_list)
        m = mean_list(kk);
        v = std_list(jj);
        v = v * v;
        a = ((1 - m) / v - 1 / m) * m * m;
        b = a * (1 / m - 1);
        if a < 1 || b < 1
            continue;
        end
            
        for ii = 1:length(N_list)
            N = N_list(ii);
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

            res(kk, jj, ii) = sum(abs(n_cdf - bb_cdf)) / (step - 1);
        end
    end
end

for ii = 1:length(N_list)
    ax(ii) = subplot(1, length(N_list), ii);
    for kk = 1:length(mean_list)
        plot(std_list, res(kk, :, ii))
        hold on;
    end
    title("N = " + num2str(N_list(ii)));
    xlim([0.05 0.3])
end

legend("0.1", "0.2", "0.3", "0.4", "0.5")
leg = legend('show');
title(leg,'Mean')

linkaxes(ax,'y');
yticklabels(ax(2:end),{})
t.TileSpacing = 'compact';
ylim([0 0.03]);