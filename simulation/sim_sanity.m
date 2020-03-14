%% global parameters
rng(0);

M1_list = 5:12;
M2_list = 5:12;

ITER = 10000;
alpha = 0.05;

%% parameters for 1
m1 = 0.52;
v1 = 0.03;

v1 = v1 * v1;
a1 = ((1 - m1) / v1 - 1 / m1) * m1 * m1;
b1 = a1 * (1 / m1 - 1);

N1 = 1000;


%% parameters for 2
m2 = 0.46;
v2 = 0.03;

v2 = v2 * v2;
a2 = ((1 - m2) / v2 - 1 / m2) * m2 * m2;
b2 = a2 * (1 / m2 - 1);

N2 = 1000;

%% Simulate 1
mat = zeros(length(M1_list), length(M2_list));
mat2 = zeros(length(M1_list), length(M2_list));
res = zeros(ITER, 1);

for kk = 1:length(M1_list)
    for jj = 1:length(M2_list)
        a = [a2 a1];
        b = [b2 b1];
        N = [N2 N1];
        M = [M2_list(jj) M1_list(kk)];
        
        Ep = a ./ (a + b);
        Vp = a .* b .* (a + b + N) ./ ((N .* (a + b) .^ 2) .* (a + b + 1));

        S = sqrt(sum((M - 1) .* Vp) ./ (sum(M) - 2));
        %nu = sum(M) - 2;
        nu = sum(Vp ./ M) ^ 2 / sum((Vp ./ M) .^ 2 ./ (M - 1));
        Et = diff(Ep) / S / sqrt(sum(1 ./ M));

        t_star = tinv(1 - alpha ./ 2, nu);

        mat2(kk, jj) = normcdf(t_star - Et); % false negative rate for 1-sided and 2-sided t-tests
    end
end

for kk = 1:length(M1_list)
    for jj = 1:length(M2_list)
        M1 = M1_list(kk);
        M2 = M2_list(jj);
        
        a = [a2 a1];
        b = [b2 b1];
        N = [N2 N1];
        M = [M2_list(jj) M1_list(kk)];
        
        Ep = a ./ (a + b);
        Vp = a .* b .* (a + b + N) ./ ((N .* (a + b) .^ 2) .* (a + b + 1));
        
        p1 = normrnd(Ep(2), sqrt(Vp(2)), M1, ITER) / N1;
        p2 = normrnd(Ep(1), sqrt(Vp(1)), M2, ITER) / N2;
        
        mat(kk, jj) = 1 - sum(P < alpha) / ITER;
    end
end

((mat2 - mat) ./ mat)'