% levy.m
% 该函数通过 Mantegna 算法生成符合 Lévy 分布的随机步长
function o = levy(d, beta)
    % d: 维度
    % beta: Lévy 指数, 通常取 1.5

    num = gamma(1 + beta) * sin(pi * beta / 2);
    den = gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2);

    sigma_u = (num / den)^(1 / beta);

    u = random('Normal', 0, sigma_u, 1, d);
    v = random('Normal', 0, 1, 1, d);

    step = u ./ (abs(v).^(1 / beta));

    o = step;
end