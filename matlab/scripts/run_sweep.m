function run_sweep(cfg)
%RUN_SWEEP  Solve every PCM x k-scenario x h batch and cache the results.
%   Each batch holds every (EG loading, source case) pair. Finished batches
%   are skipped, so an interrupted sweep can be resumed.

grid = make_grid(cfg.geom, cfg.r_in, cfg.r_out, cfg.N);
for ip = 1:numel(cfg.pcms)
    pcm = pcm_properties(cfg.pcms{ip});
    for ik = 1:numel(cfg.k_scenarios)
        ks = cfg.k_scenarios{ik};
        for h = cfg.h_list
            if h ~= cfg.h_list(1) && ~ismember(ks, cfg.h_sens_k), continue; end
            fname = fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_%s.mat', pcm.name, ks, h, cfg.mode));
            if exist(fname, 'file')
                fprintf('skip  %s\n', fname); continue
            end
            [props, src, meta] = build_batch(pcm, k_scenario_beta(ks), cfg);
            tic;
            res = solve_pcm_1d(props, grid, h, pcm.Tm - cfg.dT_sub, src, cfg.t_end, cfg.dt, cfg.record_every);
            el = toc;
            bal = max(abs(res.E(:,end) - res.Qin(:,end)) ./ abs(res.Qin(:,end)));
            res = rmfield(res, {'T', 'f'});             % keep the file small
            save(fname, 'res', 'meta', 'props', '-v7.3');
            fprintf('%-10s %-4s h=%-5d B=%d  %.0f s  iters %.2f/%d  energy balance %.1e\n', ...
                pcm.name, ks, h, numel(meta.w), el, res.mean_iter, res.max_iter, bal);
        end
    end
end

% bypass-control batches at the base heat transfer coefficient
h = cfg.h_list(1);
for r = 1:size(cfg.bypass_runs, 1)
    pcm = pcm_properties(cfg.bypass_runs{r, 1}); ks = cfg.bypass_runs{r, 2};
    fname = fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_bypass_%s.mat', pcm.name, ks, h, cfg.mode));
    if exist(fname, 'file'), fprintf('skip  %s\n', fname); continue; end
    [props, src, meta] = build_batch(pcm, k_scenario_beta(ks), cfg);
    tic;
    res = solve_pcm_1d(props, grid, h, pcm.Tm - cfg.dT_sub, src, cfg.t_end, cfg.dt, cfg.record_every, struct('bypass', true));
    el = toc;
    res = rmfield(res, {'T', 'f'});
    save(fname, 'res', 'meta', 'props', '-v7.3');
    fprintf('%-10s %-4s h=%-5d BYPASS  %.0f s\n', pcm.name, ks, h, el);
end

% natural-convection-bound batches at the base heat transfer coefficient
for r = 1:size(cfg.conv_runs, 1)
    pcm = pcm_properties(cfg.conv_runs{r, 1}); ks = cfg.conv_runs{r, 2};
    fname = fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_conv_%s.mat', pcm.name, ks, h, cfg.mode));
    if exist(fname, 'file'), fprintf('skip  %s\n', fname); continue; end
    [props, src, meta] = build_batch(pcm, k_scenario_beta(ks), cfg);
    conv = convection_params(pcm.name, meta.w, cfg.conv_w_s);
    tic;
    res = solve_pcm_1d(props, grid, h, pcm.Tm - cfg.dT_sub, src, cfg.t_end, cfg.dt, cfg.record_every, struct('conv', conv));
    el = toc;
    res = rmfield(res, {'T', 'f'});
    save(fname, 'res', 'meta', 'props', '-v7.3');
    fprintf('%-10s %-4s h=%-5d CONVECTION  %.0f s\n', pcm.name, ks, h, el);
end
end
