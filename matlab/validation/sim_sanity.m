clear;

param_sanity

ITER = 10000;

%%
mat = zeros(length(M1_list), length(M0_list));

%% Simulations
for kk = 1:length(M1_list)
    for jj = 1:length(M0_list)
        M1 = M1_list(kk);
        M0 = M0_list(jj);
        
        %a = [a0 a1];
        %b = [b0 b1];
        %N = [N0 N1];
        %M = [M0_list(jj) M1_list(kk)];
        
        p1 = betarnd(a1, b1, M1, ITER);
        p0 = betarnd(a0, b0, M0, ITER);
        
        p1 = binornd(N1, p1) / N1;
        p0 = binornd(N0, p0) / N0;
        
        s1 = std(p1);
        s0 = std(p0);
        
        %nu = (s1 .^ 2 / M1 + s0 .^ 2 / M0) .^ 2 / ((s1 .^ 2 / M1) .^ 2 / (M1 - 1) + (s0 .^ 2 / M0) .^ 2 / (M0 - 1));
        %t = (mean(p1) - mean(p0)) ./ sqrt(s1 .^ 2 / M1 + s0 .^ 2 / M0);
        %t_star = tinv(1 - alpha ./ 2, nu);
        
        %mat(kk, jj) = 1 - sum(t > t_star) / ITER;
        [H,P] = ttest2(p1, p0, 'Vartype', 'unequal', 'Tail', 'right');
        mat(kk, jj) = 1 - sum(P < alpha) / ITER;
    end
end

mat

csvwrite('sim.csv', mat)
