function cfg = project_config()
%PROJECT_CONFIG Return project-relative locations used by the analyses.
%   Edit only this file if large simulation files are stored outside the
%   repository. Environment variables take precedence over local defaults.

root = fileparts(fileparts(mfilename('fullpath')));
cfg.root = root;
cfg.simulationData = getenv('DPF_SIMULATION_DATA');
if isempty(cfg.simulationData)
    cfg.simulationData = fullfile(root, 'data', 'simulation');
end
cfg.experimentalData = fullfile(root, 'data', 'experimental');
cfg.characteristics = fullfile(root, 'utilities', 'characteristics');
cfg.output = fullfile(root, 'results');
end
