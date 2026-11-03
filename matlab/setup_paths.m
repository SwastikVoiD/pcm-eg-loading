function setup_paths()
%SETUP_PATHS  Add all project folders to the MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'scripts'));
if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
end
