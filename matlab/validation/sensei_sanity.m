clear;

param_sanity
%%
mat = zeros(length(M1_list), length(M0_list));

%% Sensei
for kk = 1:length(M1_list)
    for jj = 1:length(M0_list)
        a = [a0 a1];
        b = [b0 b1];
        N = [N0 N1];
        M = [M0_list(jj) M1_list(kk)];
        
        Ep = a ./ (a + b);
        Vp = a .* b .* (a + b + N) ./ ((N .* (a + b) .^ 2) .* (a + b + 1));

        %nu = sum(M) - 2;
        nu = sum(Vp ./ M) ^ 2 / sum((Vp ./ M) .^ 2 ./ (M - 1));
        Et = diff(Ep) / sqrt(sum(Vp ./ M));

        t_star = tinv(1 - alpha, nu);

        % mat(kk, jj) = normcdf(t_star - Et); % false negative rate for 1-sided and 2-sided t-tests
        mat(kk, jj) = tcdf(t_star - Et, nu);
    end
end

mat

csvwrite('sensei.csv', mat)
