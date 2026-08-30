function startup()
%STARTUP Add the repository's source directories to the MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'analysis')));
addpath(genpath(fullfile(root, 'config')));
addpath(genpath(fullfile(root, 'experiments')));
addpath(genpath(fullfile(root, 'utilities')));
end
