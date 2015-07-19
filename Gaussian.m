classdef Gaussian < ProbabilityDensityFunction
    properties (Access = 'protected')
        % dataset parameters
        RVNames;
        infoStr = '';
        
        % model parameters
        parameters;
        
        % inner state
        isEstimated = 0;
    end % properties
    
    methods 
        % Constructor
        function obj = Gaussian (dataTable)
            if (nargin == 0)
                obj.isEstimated = 0;
            else 
                obj.RVNames = dataTable.Properties.VariableNames;
                data        = table2array (dataTable);
                obj.estimateModel_ (data);
                obj.isEstimated = 1;
            end % if
        end % function
        
        % Marginalization
        function obj = marg (obj, RVNames)
            checkRVNames (obj.RVNames, RVNames);
            
            [~, J] = getIJ (obj.RVNames, RVNames);
            obj.parameters.mu    = obj.parameters.mu(J);
            obj.parameters.sigma = obj.parameters.sigma(J, J);
            
            obj.RVNames = RVNames;
        end % function
        
        % Conditioning
        function obj = cond (obj, RVNames, x_J)
            checkRVNames (obj.RVNames, RVNames);
            
            [I, J] = getIJ (obj.RVNames, RVNames);
            
            % x_i - mu'_i, with mu'_i = mu_i + sig_IJ * inv(sig_JJ) * (x_J - mu_J)
            obj.parameters.mu = obj.parameters.mu(I) ...
                + obj.parameters.sigma(I, J) ...
                * inv (obj.parameters.sigma(J, J)) ...
                * (x_J - obj.parameters.mu(J));
            
            obj.parameters.sigma = upperSchurComplement (obj.parameters.sigma, I, J);
            
            obj.RVNames = obj.RVNames(~ismember (obj.RVNames, RVNames));
        end % function
        
        % Display model-summary
        function disp (obj)         
            if obj.isEstimated
                disp ('mean');
                disp (array2table (obj.parameters.mu', ...
                    'VariableNames', obj.RVNames));
                disp ('covariance-matrix');
                disp (array2table (obj.parameters.sigma, ...
                    'VariableNames', obj.RVNames, ...
                    'RowNames', obj.RVNames));
                disp ('precision-matrix');
                disp (array2table (obj.parameters.prec, ...
                    'VariableNames', obj.RVNames, ...
                    'RowNames', obj.RVNames));
            else 
                warning ('No parameters has been trained till now.');
            end % if
        end % function
        
        % Densitiy value for given data
        function prob = p (obj, X)
            prob = mvnpdf (X, obj.parameters.mu', obj.parameters.sigma);
        end % function
        
        % Plot the density
        function plot (obj, data)
            [D, ~] = size (obj.parameters.mu);
            subplot1 (D-1, D-1, 'Gap', [0.01, 0.01]);
            
            p = 1;
            if (nargin < 2) % only plot density as heatmap
                error ('TODO: Implement the Gaussian heat-map.');           
                for i = 1:D ; for j = 1:D
                    if (i ~= j)
                        subplot (ind2sub ([D-1, D-1], p));
                        p = p + 1;
                    end % if 
                end ; end % for
            else            % plot density and data joint
                for i = 2:D ; for j = 1:(i-1)     
                    % plot the data-points
                    subplot1 ((i-2) * (D-1) + j);
                    C = linspecer (1);
                    scatter (data(:, j), data(:, i), 25, C, '*');
                    
                    % put labels to the axes
                    if (j == 1)
                        ylabel (obj.RVNames(i));
                    end % if
                    if (i == D)
                        xlabel (obj.RVNames(j));
                    end % if

                    % plot gaussians
                    [mu, sigma, ~] = obj.marg_const_ ({obj.RVNames{j}, obj.RVNames{i}});
                    hold on;                        
                    error_ellipse2 (sigma, mu, 0.9, 'Color', C);
                    hold off;
                    p = p + 1;
                end ; end % for
            end % if
        end % function
    end % methods
    
    methods (Access = 'protected')
        function estimateModel_ (obj, data, ~)
            obj.parameters = struct ('mu', mean (data)', 'sigma', cov (data));
            obj.parameters.prec = inv (obj.parameters.sigma);
        end % function
        
        % Marginalize without side-effects
        function [mu, sigma, prec] = marg_const_ (obj, RVNames)
            checkRVNames (obj.RVNames, RVNames);
            
            [~, J] = getIJ (obj.RVNames, RVNames);
            
            mu    = obj.parameters.mu(J);
            sigma = obj.parameters.sigma(J, J);
            prec  = obj.parameters.prec(J, J);
        end % function
    end % methods
end % class