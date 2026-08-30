function tests = test_repository_structure
%TEST_REPOSITORY_STRUCTURE Basic checks that do not require research data.
tests = functiontests(localfunctions);
end

function configurationResolves(testCase)
startup;
cfg = project_config;
verifyTrue(testCase, isfolder(cfg.root));
verifyTrue(testCase, isfolder(cfg.simulationData));
verifyTrue(testCase, isfolder(cfg.experimentalData));
end

function entryPointsExist(testCase)
startup;
names = {'analyze_dpf_parameter_sweep', ...
         'validate_phantom_measurements', ...
         'evaluate_absolute_recovery', ...
         'evaluate_differential_recovery'};
for k = 1:numel(names)
    verifyNotEmpty(testCase, which(names{k}));
end
end
