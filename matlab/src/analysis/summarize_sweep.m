function S = summarize_sweep(res, meta, cfg, t_c)
%SUMMARIZE_SWEEP  Optimum loading and design penalty for every source type.
%   S = summarize_sweep(res, meta, cfg, t_c) returns a table, one row per
%   source (constant, then each shape/amplitude/period), with
%     w_opt       optimum EG fraction of the phase-averaged E(w) curve
%     w_opt_fs    same, constrained to w >= cfg.w_formstable
%     w_opt_lo/hi range of per-start-phase optima
%     E_opt       stored energy at the optimum [kJ/m]
%     fl_opt      melt fraction at the optimum
%     dw          w_opt - w_opt(constant)
%     penalty_pct energy lost when the unit uses w_opt(constant) instead of
%                 its own optimum:  100*(1 - E(w_opt_const)/E_opt)
%     dE_pct      change of the achievable optimum vs the constant source

[E, fl] = energy_at_window(res, t_c);
cases = meta.cases;
keys = unique(cases(:, 1:3), 'rows', 'stable');
w = cfg.w(:);
nG = size(keys, 1);
shape = strings(nG, 1);
[A, P, w_opt, w_opt_fs, w_lo, w_hi, E_opt, fl_opt, dw, penalty, dE] = deal(zeros(nG, 1));
Ebar_all = zeros(numel(w), nG);

for g = 1:nG
    cidx = find(ismember(cases(:, 1:3), keys(g, :), 'rows'));
    Ew = zeros(numel(w), numel(cidx)); Fw = Ew;
    for q = 1:numel(cidx)
        sel = meta.case == cidx(q);
        [ws, order] = sort(meta.w(sel));
        e = E(sel); fq = fl(sel);
        Ew(:, q) = e(order); Fw(:, q) = fq(order);
        [wq, ~] = optimum_loading(ws, Ew(:, q));
        if q == 1, w_lo(g) = wq; w_hi(g) = wq; end
        w_lo(g) = min(w_lo(g), wq); w_hi(g) = max(w_hi(g), wq);
    end
    Ebar = mean(Ew, 2); Fbar = mean(Fw, 2);
    Ebar_all(:, g) = Ebar;
    [w_opt(g), E_opt(g)] = optimum_loading(w, Ebar);
    w_opt_fs(g) = optimum_loading(w, Ebar, cfg.w_formstable);
    fl_opt(g) = interp1(w, Fbar, w_opt(g));
    shape(g) = shape_name(keys(g, 1));
    A(g) = keys(g, 2); P(g) = keys(g, 3);
end
% metrics relative to the constant source (row 1)
for g = 1:nG
    dw(g) = w_opt(g) - w_opt(1);
    E_at_const_opt = interp1(w, Ebar_all(:, g), w_opt(1), 'pchip');
    penalty(g) = 100 * (1 - E_at_const_opt / E_opt(g));
    dE(g) = 100 * (E_opt(g) / E_opt(1) - 1);
end
S = table(shape, A, P/60, w_opt, w_opt_fs, w_lo, w_hi, E_opt/1e3, fl_opt, dw, penalty, dE, ...
    'VariableNames', {'shape', 'A_K', 'P_min', 'w_opt', 'w_opt_fs', 'w_opt_lo', 'w_opt_hi', ...
                      'E_opt_kJ_per_m', 'fl_opt', 'dw', 'penalty_pct', 'dE_pct'});
end

function s = shape_name(id)
names = ["constant", "sine", "square"];
s = names(id + 1);
end
