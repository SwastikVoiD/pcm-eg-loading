function ALL = run_analysis(cfg)
%RUN_ANALYSIS  Optimum loading and design penalty for every batch and window.
%   Writes results/summary_all.csv (one row per PCM, k scenario, h, window, source).

ALL = table();
for ip = 1:numel(cfg.pcms)
    for ik = 1:numel(cfg.k_scenarios)
        for h = cfg.h_list
            fname = fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_%s.mat', ...
                cfg.pcms{ip}, cfg.k_scenarios{ik}, h, cfg.mode));
            if ~exist(fname, 'file'), continue; end
            L = load(fname);
            for t_c = cfg.t_windows
                S = summarize_sweep(L.res, L.meta, cfg, t_c);
                n = height(S);
                S = [table(repmat(string(cfg.pcms{ip}), n, 1), repmat(string(cfg.k_scenarios{ik}), n, 1), ...
                           repmat(h, n, 1), repmat(t_c/60, n, 1), ...
                           'VariableNames', {'pcm', 'k_scenario', 'h', 't_c_min'}), S]; %#ok<AGROW>
                ALL = [ALL; S]; %#ok<AGROW>
            end
        end
    end
end
ALL.bypass = false(height(ALL), 1);

% bypass-control batches
h = cfg.h_list(1);
for r = 1:size(cfg.bypass_runs, 1)
    fname = fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_bypass_%s.mat', ...
        cfg.bypass_runs{r,1}, cfg.bypass_runs{r,2}, h, cfg.mode));
    if ~exist(fname, 'file'), continue; end
    L = load(fname);
    for t_c = cfg.t_windows
        S = summarize_sweep(L.res, L.meta, cfg, t_c);
        n = height(S);
        S = [table(repmat(string(cfg.bypass_runs{r,1}), n, 1), repmat(string(cfg.bypass_runs{r,2}), n, 1), ...
                   repmat(h, n, 1), repmat(t_c/60, n, 1), ...
                   'VariableNames', {'pcm', 'k_scenario', 'h', 't_c_min'}), S]; %#ok<AGROW>
        S.bypass = true(n, 1);
        ALL = [ALL; S]; %#ok<AGROW>
    end
end
ALL.convection = false(height(ALL), 1);

% natural-convection-bound batches
for r = 1:size(cfg.conv_runs, 1)
    fname = fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_conv_%s.mat', ...
        cfg.conv_runs{r,1}, cfg.conv_runs{r,2}, h, cfg.mode));
    if ~exist(fname, 'file'), continue; end
    L = load(fname);
    for t_c = cfg.t_windows
        S = summarize_sweep(L.res, L.meta, cfg, t_c);
        n = height(S);
        S = [table(repmat(string(cfg.conv_runs{r,1}), n, 1), repmat(string(cfg.conv_runs{r,2}), n, 1), ...
                   repmat(h, n, 1), repmat(t_c/60, n, 1), ...
                   'VariableNames', {'pcm', 'k_scenario', 'h', 't_c_min'}), S]; %#ok<AGROW>
        S.bypass = false(n, 1);
        S.convection = true(n, 1);
        ALL = [ALL; S]; %#ok<AGROW>
    end
end
writetable(ALL, fullfile(cfg.results_dir, 'summary_all.csv'));
fprintf('%d summary rows written.\n', height(ALL));
end
