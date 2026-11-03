function cfg = default_config(mode)
%DEFAULT_CONFIG  All study parameters in one place.
%   cfg = default_config()         full study (as reported in the paper)
%   cfg = default_config('quick')  reduced sweep for testing (~20x faster)
%
%   Every number used by the scripts lives here, so a change to the study
%   design is made in exactly one file.

if nargin < 1, mode = 'full'; end

% ---------------- geometry (shell-and-tube cross-section, Li et al. 2020)
cfg.geom   = 'radial';      % 'radial' annulus or 'planar' slab
cfg.r_in   = 0.0125;        % inner (tube) radius [m]
cfg.r_out  = 0.025;         % outer (shell) radius [m], outer wall adiabatic

% ---------------- numerics (see run_independence.m)
cfg.N            = 80;      % finite volumes across the PCM layer
cfg.dt           = 2.0;     % time step [s]
cfg.t_end        = 7200;    % simulated charging horizon [s]
cfg.record_every = 60;      % output interval [s]

% ---------------- materials and composite design space
cfg.pcms        = {'erythritol', 'paraffin'};
cfg.k_scenarios = {'low', 'mid', 'high'};   % see k_scenario_beta.m
cfg.w           = 0:0.005:0.25;             % EG mass fraction sweep [-]
cfg.w_formstable = 0.10;                    % form-stability threshold [-]

% ---------------- heat source
cfg.dT_super = 20;          % mean HTF temperature above melting point [K]
cfg.dT_sub   = 10;          % initial PCM temperature below melting point [K]
cfg.h_list   = [1000 200 5000];            % HTF-side coefficients [W/m2K]; first = base
cfg.h_sens_k = {'low'};                    % k routes that get the h sensitivity (others: base h only)
cfg.periods  = [600 1800 3600];            % source periods [s]
cfg.amps     = [10 20];                    % source amplitudes [K]
cfg.phases   = [0 pi/2 pi 3*pi/2];         % start phases averaged out [rad]

% ---------------- bypass control study (HTF stops when colder than the store)
cfg.bypass_runs = {'erythritol', 'low'; 'erythritol', 'mid'; 'erythritol', 'high'; 'paraffin', 'mid'};

% ---------------- natural-convection bound (upper bound; conduction-only = lower bound)
cfg.conv_runs = {'erythritol', 'low'; 'erythritol', 'mid'; 'erythritol', 'high'; ...
                 'paraffin', 'low'; 'paraffin', 'mid'; 'paraffin', 'high'};
cfg.conv_w_s  = 0.10;     % EG loading that fully suppresses melt convection

% ---------------- design objective
cfg.t_windows = [1800 3600 7200];          % charging windows [s]; 3600 = base

% ---------------- folders
root = fileparts(fileparts(mfilename('fullpath')));
cfg.root        = root;
cfg.results_dir = fullfile(root, 'results');
cfg.figures_dir = fullfile(root, 'figures');

if strcmpi(mode, 'quick')
    cfg.w       = 0:0.01:0.25;
    cfg.h_list  = 1000;
    cfg.phases  = [0 pi];
    cfg.k_scenarios = {'mid'};
    cfg.bypass_runs = {'erythritol', 'mid'};
    cfg.conv_runs   = {'erythritol', 'mid'};
end
cfg.mode = mode;
end
