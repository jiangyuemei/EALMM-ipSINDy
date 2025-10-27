%% =========================================================================
% EABDF-ipSINDy: Ensemble Adaptive BDF with Inner-Product Sparse Identification
% Robust Nonlinear Dynamics Discovery from Noisy Observational Data
% 
% Computational Framework for:
% - Multi-order Backward Differentiation Formulae (BDF) Methods
% - Ensemble-based Adaptive Denoising
% - Inner-Product Driven Sparse Regression
% - Cross-validated Model Selection
%
% Target System: Lorenz Attractor (Chaotic Regime)
% =========================================================================

clear all; close all; clc;

%% =========================================================================
% EXPERIMENTAL CONFIGURATION
% =========================================================================

% Noise-to-Signal Ratio Spectrum
noiseLevels = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1]; 

% BDF Method Orders for Comparative Analysis
bdfOrders = 1:5;  % BDF-1 through BDF-5

% Performance Metric Storage
errorMatrix = zeros(length(bdfOrders), length(noiseLevels));

%% =========================================================================
% REFERENCE SYSTEM: LORENZ ATTRACTOR PARAMETERIZATION
% =========================================================================

systemParameters = [10; 28; 8/3]; % Chaotic regime parameters
systemDimension = 3;
initialState = [-8; 7; 27];  % Canonical initial conditions

% Temporal Discretization
timeStep = 0.01;
timeDomain = 0:timeStep:10;
integrationOptions = odeset('RelTol', 1e-10, 'AbsTol', 1e-10 * ones(1, 3));

% Generate Reference Trajectories
[timeVector, stateTrajectories] = ode45(@(t, x) lorenz(t, x, systemParameters), ...
                                       timeDomain, initialState, integrationOptions);
[dataPoints, stateDimensions] = size(stateTrajectories);

%% =========================================================================
% GROUND TRUTH SPECIFICATION
% =========================================================================

libraryCardinality = 10;
groundTruthCoefficients = zeros(libraryCardinality, 3);

% Lorenz System Coefficient Mapping
groundTruthCoefficients(2, 1) = -10;   % \dot{u}_1: -10u_1
groundTruthCoefficients(2, 2) = 28;    % \dot{u}_2: 28u_1  
groundTruthCoefficients(3, 1) = 10;    % \dot{u}_1: 10u_2
groundTruthCoefficients(3, 2) = -1;    % \dot{u}_2: -u_2
groundTruthCoefficients(4, 3) = -8/3;  % \dot{u}_3: -8/3 u_3
groundTruthCoefficients(6, 3) = 1;     % \dot{u}_3: u_1u_2
groundTruthCoefficients(7, 2) = -1;    % \dot{u}_2: -u_1u_3

%% =========================================================================
% COMPUTATIONAL FRAMEWORK EXECUTION
% =========================================================================

fprintf('Initiating EABDF-ipSINDy Analysis Framework...\n\n');

for orderIndex = 1:length(bdfOrders)
    currentOrder = bdfOrders(orderIndex);
    fprintf('Processing BDF-%d Method Configuration...\n', currentOrder);
    
    for noiseIndex = 1:length(noiseLevels)
        currentNoiseLevel = noiseLevels(noiseIndex);
        fprintf('  Noise Level: σ_NR = %.0e\n', currentNoiseLevel);
        
        % -----------------------------------------------------------------
        % BDF OPERATOR CONSTRUCTION
        % -----------------------------------------------------------------
        [A_matrix, B_matrix] = BDF(length(timeVector) - 1, currentOrder);
        
        % -----------------------------------------------------------------
        % ENSEMBLE CONFIGURATION
        % -----------------------------------------------------------------
        ensembleReplicates = 10;
        polynomialOrder = 2;
        librarySize = 10;
        coefficientEnsemble = zeros(librarySize, systemDimension, ensembleReplicates);
        stateEnsemble = zeros(length(timeVector), systemDimension, ensembleReplicates);
        
        %% ENSEMBLE PROCESSING PIPELINE
        for ensembleIndex = 1:ensembleReplicates
            % Noise Injection
            noiseSigma = currentNoiseLevel * norm(stateTrajectories, 'fro') / ...
                        sqrt(dataPoints * stateDimensions);
            noisyStates = stateTrajectories + normrnd(0, noiseSigma, size(stateTrajectories));
            
            % -----------------------------------------------------------------
            % ADAPTIVE MOVING AVERAGE FILTERING
            % -----------------------------------------------------------------
            smoothingWindow = ceil(length(timeDomain) / 100);  
            for smoothingIteration = 1:200
                filterKernels = cell(systemDimension, 1); 
                for stateIndex = 1:systemDimension
                    [filterKernels{stateIndex}, ~, ~] = ...
                        AMAF(timeDomain, noisyStates(:, stateIndex), [], ...
                        smoothingWindow, [], [], [], [], []);
                    
                    kernelRadius = (length(filterKernels{stateIndex}) - 1) / 2;
                    symmetricExtension = [flipud(noisyStates(2:kernelRadius+1, stateIndex)); 
                                        noisyStates(:, stateIndex);
                                        flipud(noisyStates(end-kernelRadius:end-1, stateIndex))];
                    stateEnsemble(:, stateIndex, ensembleIndex) = ...
                        conv(symmetricExtension, filterKernels{stateIndex}, 'valid');
                end
            end    
        end
        
        % Ensemble Averaging
        denoisedStates = mean(stateEnsemble, 3);
        
        % -----------------------------------------------------------------
        % FEATURE LIBRARY CONSTRUCTION
        % -----------------------------------------------------------------
        featureLibrary = poolData(denoisedStates, systemDimension, polynomialOrder);
        
        %% SPARSE IDENTIFICATION - FIRST STAGE
        sparsityParameter = 0.25;
        identifiedCoefficients = sparsifyDynamics(timeStep * B_matrix * featureLibrary, ...
                                                A_matrix * denoisedStates, ...
                                                sparsityParameter, systemDimension);
        
        coefficientEnsemble(:, :, ensembleIndex) = identifiedCoefficients;
        
        % Ensemble Coefficient Aggregation
        aggregatedCoefficients = mean(coefficientEnsemble, 3);
        
        %% INNER-PRODUCT-DRIVEN REFINEMENT
        refinedCoefficients = zeros(librarySize, systemDimension);
        
        % Non-zero Coefficient Identification
        nonzeroIndices = cell(systemDimension, 1);
        for dim = 1:systemDimension
            nonzeroIndices{dim} = find(aggregatedCoefficients(:, dim));
        end
        
        % Cross-Validation Framework
        crossValidationParameter = 0.25;
        trainingRatio = 0.8;
        trainingSamples = floor(trainingRatio * size(featureLibrary, 1));
        trainingIndices = 1:trainingSamples;
        validationIndices = trainingSamples + 1:size(featureLibrary, 1);
        
        % BDF Operators for Cross-Validation
        trainingSteps = length(trainingIndices) - 1; 
        [A_train, B_train] = BDF(trainingSteps, currentOrder);
        validationSteps = length(validationIndices) - 1; 
        [A_val, B_val] = BDF(validationSteps, currentOrder);
        
        % Optimal Sparsity Level Selection
        minimumError = inf(1, systemDimension);
        
        for dimensionIndex = 1:systemDimension
            for sparsityLevel = 1:length(nonzeroIndices{dimensionIndex})
                % Sparse Regression with Fixed Sparsity
                sparseCoefficients = ipSINDy_K(...
                    timeStep * B_train * featureLibrary(trainingIndices, nonzeroIndices{dimensionIndex}), ...
                    A_train * denoisedStates(trainingIndices, dimensionIndex), ...
                    crossValidationParameter, sparsityLevel, 1);
                
                % Validation Residual Computation
                validationResidual = timeStep * B_val * ...
                    featureLibrary(validationIndices, nonzeroIndices{dimensionIndex}) * sparseCoefficients - ...
                    A_val * denoisedStates(validationIndices, dimensionIndex);
                residualNorm = norm(validationResidual);
                
                % Model Selection Criterion
                if residualNorm < minimumError(dimensionIndex)
                    minimumError(dimensionIndex) = residualNorm;
                    refinedCoefficients(:, dimensionIndex) = 0;
                    refinedCoefficients(nonzeroIndices{dimensionIndex}, dimensionIndex) = sparseCoefficients;
                end
            end
        end
        
        % Performance Quantification
        errorMatrix(orderIndex, noiseIndex) = ...
            norm(refinedCoefficients - groundTruthCoefficients, 'fro') / ...
            norm(groundTruthCoefficients, 'fro');
    end
end

%% =========================================================================
% COMPUTATIONAL RESULTS
% =========================================================================

fprintf('\n=== EABDF-ipSINDy IDENTIFICATION PERFORMANCE ===\n');
disp('Relative Coefficient Errors (Normalized Frobenius Norm):');
disp(errorMatrix);

%% =========================================================================
% PERFORMANCE VISUALIZATION
% =========================================================================

figure('Position', [100, 100, 1200, 800]); hold on;

markerStyles = {'-o', '-s', '-^', '-d', '-p'};
colorPalette = lines(length(bdfOrders));

for methodIndex = 1:length(bdfOrders)
    loglog(noiseLevels, errorMatrix(methodIndex, :), markerStyles{methodIndex}, ...
           'LineWidth', 2.2, 'MarkerSize', 8, 'Color', colorPalette(methodIndex, :));
end

% Visualization Enhancement
legendConfiguration = {
    'EABDF-1-ipSINDy',...
    'EABDF-2-ipSINDy',...
    'EABDF-3-ipSINDy',...
    'EABDF-4-ipSINDy',...
    'EABDF-5-ipSINDy'
};

legend(legendConfiguration, 'Location', 'southwest', 'FontSize', 11);
xlabel('\sigma_{NR}', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('e(C)', 'FontSize', 13, 'FontWeight', 'bold');
% title('EABDF-ipSINDy: Robust Identification Performance', 'FontSize', 14, 'FontWeight', 'bold');

grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
box on;

fprintf('\nAnalysis completed successfully.\n');