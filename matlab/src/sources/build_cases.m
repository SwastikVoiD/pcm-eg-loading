function cases = build_cases(cfg)
%BUILD_CASES  Table of source cases: one row per [shape_id, A, P, phase].
%   Row 1 is the constant source; then sine and square waves for every
%   period, amplitude and start phase in cfg.

cases = [0 0 3600 0];
for sid = [1 2]
    for P = cfg.periods
        for A = cfg.amps
            for ph = cfg.phases
                cases(end+1, :) = [sid A P ph]; %#ok<AGROW>
            end
        end
    end
end
end
