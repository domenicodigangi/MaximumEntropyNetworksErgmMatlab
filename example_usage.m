%% Maximum Entropy Network Reconstruction - Example Usage
% This script demonstrates how to use the Max_Entr_Nets package for
% network reconstruction and systemic risk analysis.
%
% Author: Domenico Di Gangi
% Date: 2018

%% Example 1: Basic Network Reconstruction with BIPWCM
% Bipartite Weighted Configuration Model - estimates network from strength sequences

fprintf('\n=== Example 1: Bipartite Weighted Configuration Model ===\n\n');

% Define a simple bipartite network (e.g., 4 investors, 3 assets)
% Rows represent investors, columns represent assets
true_network = [50, 30, 20;
                40, 60, 50;
                80, 70, 50;
                30, 40, 50];

% Extract strength sequences (row and column sums)
investor_strengths = sum(true_network, 2);  % Row sums
asset_strengths = sum(true_network, 1)';    % Column sums

fprintf('Original network:\n');
disp(true_network);
fprintf('Investor strengths: %s\n', mat2str(investor_strengths'));
fprintf('Asset strengths: %s\n\n', mat2str(asset_strengths'));

% Prepare input data for the model
in_data = {investor_strengths, asset_strengths};

% Estimate the Bipartite Weighted Configuration Model
fprintf('Estimating BIPWCM model...\n');
model_BIPWCM = Max_Entr_Nets('BIPWCM', in_data);

fprintf('Model estimated in %.4f seconds\n', model_BIPWCM.estimation_time);
fprintf('Maximum relative error: %.6f\n\n', max(model_BIPWCM.precision.errors));

% Get the expected adjacency matrix
expected_matrix = model_BIPWCM.exp_matrix();
fprintf('Expected adjacency matrix:\n');
disp(expected_matrix);

% Sample networks from the ensemble
n_samples = 10;
fprintf('Sampling %d networks from the ensemble...\n', n_samples);
network_samples = model_BIPWCM.sample(n_samples);

fprintf('Sample 1:\n');
disp(network_samples{1});


%% Example 2: Enhanced Configuration Model (BIPECM)
% Uses both strength and degree sequences for better reconstruction

fprintf('\n=== Example 2: Enhanced Configuration Model ===\n\n');

% Binary version of the network (for degree sequences)
binary_network = (true_network > 0);
investor_degrees = sum(binary_network, 2);
asset_degrees = sum(binary_network, 1)';

% Prepare input data with both strength and degree sequences
in_data_enhanced = {investor_strengths, asset_strengths, ...
                    investor_degrees, asset_degrees};

fprintf('Estimating BIPECM model...\n');
model_BIPECM = Max_Entr_Nets('BIPECM', in_data_enhanced);

fprintf('Model estimated in %.4f seconds\n', model_BIPECM.estimation_time);
fprintf('Maximum relative error: %.6f\n\n', max(max(model_BIPECM.precision.errors)));

expected_matrix_ecm = model_BIPECM.exp_matrix();
fprintf('Expected adjacency matrix (BIPECM):\n');
disp(expected_matrix_ecm);


%% Example 3: MECAPM Model
% Maximum Entropy CAPM - based on capitalization

fprintf('\n=== Example 3: Maximum Entropy CAPM Model ===\n\n');

% For MECAPM, we use the strength sequences
in_data_capm = {investor_strengths, asset_strengths};

fprintf('Estimating MECAPM model...\n');
model_MECAPM = Max_Entr_Nets('MECAPM', in_data_capm);

fprintf('Model estimated in %.4f seconds\n', model_MECAPM.estimation_time);

expected_matrix_capm = model_MECAPM.exp_matrix();
fprintf('Expected adjacency matrix (MECAPM):\n');
disp(expected_matrix_capm);


%% Example 4: Systemic Risk Analysis - Vulnerable Banks Framework
% Demonstrates fire sales spillover analysis

fprintf('\n=== Example 4: Systemic Risk Analysis ===\n\n');

% Define equity values for each investor
equity_values = [20; 30; 40; 25];  % Equity for each investor

% Define shock scenario (1% depreciation on all assets)
shock_scenario = [0.01, 0.01, 0.01];

fprintf('Computing systemic risk with REAL network...\n');
[AV_real, SYS_real, VUL_real] = Vulnerable_Banks('REAL', ...
                                                  true_network, ...
                                                  equity_values, ...
                                                  shock_scenario);

fprintf('Aggregate Vulnerability: %.4f\n', AV_real);
fprintf('Systemicness by investor:\n');
disp(SYS_real);
fprintf('Vulnerability by investor:\n');
disp(VUL_real);

% Now estimate risk from partial information using BIPWCM
fprintf('\nEstimating systemic risk from PARTIAL information (BIPWCM)...\n');
[AV_est, SYS_est, VUL_est] = Vulnerable_Banks('ESTIMATE-BIPWCM', ...
                                               in_data, ...
                                               equity_values, ...
                                               shock_scenario);

fprintf('Estimated Aggregate Vulnerability: %.4f\n', AV_est);
fprintf('Estimated Systemicness by investor:\n');
disp(SYS_est);

% Compare results
fprintf('\n--- Comparison ---\n');
fprintf('AV (Real): %.4f vs AV (Estimated): %.4f\n', AV_real, AV_est);
fprintf('Relative error: %.2f%%\n', 100*abs(AV_real - AV_est)/AV_real);


%% Example 5: Listing Available Models

fprintf('\n=== Example 5: Available Models ===\n\n');

fprintf('All available models:\n');
Max_Entr_Nets('LIST');

fprintf('\nBipartite models only:\n');
Max_Entr_Nets('LIST-BIP');


%% Example 6: Precision Control

fprintf('\n=== Example 6: Precision Control ===\n\n');

% Set higher precision requirement
high_precision = 1e-4;
fprintf('Estimating BIPWCM with precision = %.0e\n', high_precision);

model_precise = Max_Entr_Nets('BIPWCM', in_data, high_precision);
fprintf('Achieved maximum error: %.6f\n', max(model_precise.precision.errors));
fprintf('Mean error: %.6f\n', mean(model_precise.precision.errors));


%% Summary

fprintf('\n========================================\n');
fprintf('Examples completed successfully!\n');
fprintf('========================================\n\n');

fprintf('Key takeaways:\n');
fprintf('1. BIPWCM: Reconstruct weighted networks from strength sequences\n');
fprintf('2. BIPECM: Use both strength and degree for better accuracy\n');
fprintf('3. MECAPM: CAPM-based reconstruction from capitalization\n');
fprintf('4. Vulnerable Banks: Assess systemic risk with complete or partial network data\n');
fprintf('5. Precision control: Adjust optimization tolerance as needed\n\n');

fprintf('For more details, see:\n');
fprintf('- Individual model files in models/ directory\n');
fprintf('- Vulnerable_Banks.m for systemic risk options\n');
fprintf('- readme.pdf for technical documentation\n\n');
