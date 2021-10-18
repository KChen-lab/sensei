%% global parameters
rng(0);

M1_list = 5:12;
M2_list = 5:12;

ITER = 1000;
alpha = 0.05;

%% parameters for 1
m1 = 0.5;
v1 = 0.1;

v1 = v1 * v1;
a1 = ((1 - m1) / v1 - 1 / m1) * m1 * m1;
b1 = a1 * (1 / m1 - 1);

N1 = 1000;


%% parameters for 2
m2 = 0.3;
v2 = 0.15;

v2 = v2 * v2;
a2 = ((1 - m2) / v2 - 1 / m2) * m2 * m2;
b2 = a2 * (1 / m2 - 1);

N2 = 1000;

%% Simulate 1
mat = zeros(length(M1_list), length(M2_list));
mat2 = zeros(length(M1_list), length(M2_list));
res = zeros(ITER, 1);


%% Prediction
for kk = 1:length(M1_list)
    for jj = 1:length(M2_list)
        a = [a2 a1];
        b = [b2 b1];
        N = [N2 N1];
        M = [M2_list(jj) M1_list(kk)];
        
        Ep = a ./ (a + b);
        Vp = a .* b .* (a + b + N) ./ ((N .* (a + b) .^ 2) .* (a + b + 1));

        nu = sum(Vp ./ M) ^ 2 / sum((Vp ./ M) .^ 2 ./ (M - 1));
        Et = diff(Ep) / sqrt(sum(Vp ./ M));

        t_star = tinv(1 - alpha ./ 2, nu);

        mat2(kk, jj) = normcdf(t_star - Et); % false negative rate for 1-sided and 2-sided t-tests
    end
end


%% Simulation
for kk = 1:length(M1_list)
    for jj = 1:length(M2_list)
        M1 = M1_list(kk);
        M2 = M2_list(jj);
        p1 = binornd(N1, betarnd(a1, b1, M1, ITER)) / N1;
        p2 = binornd(N2, betarnd(a2, b2, M2, ITER)) / N2;
        
        m1 = mean(p1);
        m2 = mean(p2);
        s1 = std(p1) .^ 2;
        s2 = std(p2) .^ 2;
        
        nu = (s2 / M2 + s1 / M1) .^ 2 ./ ...
            ((s2 / M2) .^ 2 / (M2 - 1) + (s1 / M1) .^ 2 / (M1 - 1));
        
        T = (m1 - m2) ./ sqrt(s1 / M1 + s2 / M2);
        P = tcdf(T, nu, 'upper');
        mat(kk, jj) = 1 - sum(P < alpha / 2) / ITER;
    end
end

((mat2 - mat) ./ mat)'