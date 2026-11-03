function ALL = load_summary(cfg)
%LOAD_SUMMARY  Read results/summary_all.csv with correctly typed columns.
%   ALL = load_summary(cfg)
%   writetable stores the logical columns "bypass" and "convection" as 1/0, so readtable returns
%   it as double (or as text, depending on the MATLAB release). This converts
%   it back to logical whichever way it arrives.

ALL = readtable(fullfile(cfg.results_dir, 'summary_all.csv'), 'TextType', 'string');
for name = ["bypass", "convection"]
    if ~ismember(name, string(ALL.Properties.VariableNames))
        ALL.(name) = false(height(ALL), 1); continue
    end
    b = ALL.(name);
    if isnumeric(b) || islogical(b)
        b = logical(b);
    else
        b = ismember(lower(string(b)), ["true", "1"]);
    end
    ALL.(name) = b;
end
end
