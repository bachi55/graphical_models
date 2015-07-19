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
        % Input: 
        %   - dataTable ... matrix of type Datatable containing the RV 
        %                   names as column-names
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
        
        % Returns the probability for given data
        % Input:
        %   - X ... data-matrix containing data-points row-wise
        function prob = p (obj, X)
            prob = mvnpdf (X, obj.parameters.mu', obj.parameters.sigma);
        end % function
        
        % Plot the density
        % Calls:
        %   - plot():               plots the densities with error-ellipses
        %   - plot(data):           plots the densities joint with the data
        %   - plot(data, labels):   colors the data-points according to
        %                           there label
        % Input:
        %   - data   ... matrix of type Datatable
        %   - labels ... (numerical) vector containing the labels for 
        %                data-points
        function plot (obj, varargin)
            [D, K] = size (obj.parameters.mu);
            if (D > 2)
                subplot1 (D-1, D-1, 'Gap', [0.01, 0.01]);
            end % if 
            
            p = 1;
            % plot density and data joint
            
            if (nargin > 1)
                data = table2array (varargin{1}(:, obj.RVNames));
            end % if

            for i = 2:D ; for j = 1:(i-1)
                % plot the data-points
                if (D > 2)
                    subplot1 ((i-2) * (D-1) + j);
                end % if
                C = linspecer (K); % get K different colors                    
                    
      
                if (nargin > 2)     % labels are given
                    C = linspecer (numel (unique (varargin{2}))); % get K different colors
                    pointColors = C(labels, :);
                elseif (nargin > 1) % data-points are given
                    % determine the most-likely mixture component for every
                    % data-point
                    [~, mostLikelyComponent] = max (obj.p_perComponent_ (data), [], 2);    
                    pointColors = C(mostLikelyComponent, :);
                end % if

                if (nargin > 1)
                    scatter (data(:, j), data(:, i), 25, pointColors, '.');
                    hold on;
                end % if

                % put labels to the axes
                if (j == 1)
                    ylabel (obj.RVNames(i));
                end % if
                if (i == D)
                    xlabel (obj.RVNames(j));
                end % if

                % plot gaussians                       
                copiedObj = obj.marg_const ( ...
                    {obj.RVNames{j}, obj.RVNames{i}});
                [~, mu, sigma, ~] = copiedObj.getDensityParameters();

                for k = 1:K
                    plot (mu(1, k), mu(2, k), '*', 'Color', C(k, :));
                    hold on;
                    error_ellipse2 (sigma(:, :, k), mu(:, k), 0.9, ...
                        'Color', C(k, :));    
                end % for
                hold off;
                p = p + 1;
            end ; end % for
        end % function
    end % methods
    
    methods (Access = 'protected')
        function estimateModel_ (obj, data, ~)
            obj.parameters = struct ('mu', mean (data)', 'sigma', cov (data));
            obj.parameters.prec = inv (obj.parameters.sigma);
        end % function
        
        % Marginalize without side-effects
        function [mu, sigma, prec] = marg_const (obj, RVNames)
            checkRVNames (obj.RVNames, RVNames);
            
            [~, J] = getIJ (obj.RVNames, RVNames);
            
            mu    = obj.parameters.mu(J);
            sigma = obj.parameters.sigma(J, J);
            prec  = obj.parameters.prec(J, J);
        end % function
    end % methods
end % class