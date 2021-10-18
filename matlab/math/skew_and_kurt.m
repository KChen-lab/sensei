m_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];
cv_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];

res = zeros(length(m_list), length(cv_list));

for i = 1:length(m_list)
    for j = 1:length(cv_list)
        m = m_list(i);
        cv = cv_list(j);
        v = m * cv;
        v = v .* v;
        a = ((1 - m) ./ v - 1 ./ m) * m .^ 2;
        b = a * (1 / m - 1);
        if a < 1
            error('a < 1');
        end
        if b < 1 
            error('b < 1');
        end
        res(i, j) = bb_skew(10000, a, b);
    end
end
        
surf(cv_list, m_list, res)
colorbar
ylabel('mu')
xlabel('cv')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

function res = bb_skew(n, a, b)
    v = a + b;
    res = (v + 2 * n) .* (b - a) ./ (v + 2) .* sqrt((1 + v) / (n * a * b) * (n + v));
end

function res = bb_kurt(n, a, b)
    v = a + b;
    res = v .^ 2 * (1 + v) / (n .* a .* b .* (v + 2) .* (v + 3) .* (v + n));
    res = res .* (v .* (v - 1 + 6 * n) + 3 * a .* b .* (n - 2) + 6 * n .^ 2 - 3 * a .* b .* n .* (6 - n) ./ v - 18 .* a .* b .* n ^ 2 / v .^ 2);
end