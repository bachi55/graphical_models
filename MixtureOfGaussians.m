classdef MixtureOfGaussians < ProbabilityDensityFunction
	properties (Access = 'protected')
        % dataset parameters
        RVNames;
        infoStr = '';
        
        % model parameters
        parameters;
        
        % inner state
        isEstimated = 0;
        
        % debug parameters
        verbose;
    end % properties

    % Constructor
    methods 
        function obj = MixtureOfGaussians (varargin)
            % default-constructor
            if (nargin == 0)
                obj.isEstimated = 0;
            % constructor which estimates the density using EM based on the
            % data
            elseif (istable (varargin{1}))
                if (nargin < 2)
                    optPar = struct (    ...
                        'components', 1, ... % amount of mixture-components
                        'maxIt', 500,    ... % max. iterations for the EM
                        'lhTol', 1e-6,   ... % threshold for the relative lh-change (EM)
                        'verbose', 1);       % output-level: 0=no, 1=info, 2=deb
                else
                    optPar = varargin{2};
                end % if
                
                obj.verbose     = optPar.verbose;
                obj.RVNames     = varargin{1}.Properties.VariableNames;
                data            = table2array (varargin{1});
                obj.estimateModel_ (data, optPar);
                obj.isEstimated = 1;
            % copy-constructor
            elseif (isa (varargin{1}, 'ProbabilityDensityFunction'))
                rhs = varargin{1};
                
                obj.RVNames     = rhs.RVNames;
                obj.infoStr     = rhs.infoStr;
                obj.parameters  = rhs.parameters;
                obj.isEstimated = rhs.isEstimated;
                obj.verbose     = rhs.verbose;
            end % if
        end % function
        
        % Marginalization 
        % Input: 
        %   - RVNames ... contains the RVs which should be kept
        function obj = marg (obj, RVNames)
            checkRVNames (obj.RVNames, RVNames);
            
            [~, J] = getIJ (obj.RVNames, RVNames);
            
            obj.parameters.alpha = obj.parameters.alpha;
            obj.parameters.mu    = obj.parameters.mu(J, :);
            obj.parameters.sigma = obj.parameters.sigma(J, J, :);
            obj.parameters.prec  = obj.parameters.prec(J, J, :);
            
            obj.RVNames = RVNames;
        end % function
        
        % Conditioning
        function obj = cond (obj, RVNames, x_J)
            copiedObj = MixtureOfGaussians (obj);
            
            checkRVNames (copiedObj.RVNames, RVNames);
            
            [I, J] = getIJ (copiedObj.RVNames, RVNames);
            
            % get new alpha
            copiedObj.marg (RVNames);
            alpha = copiedObj.parameters.alpha .* copiedObj.p_perComponent_ (x_J')' ...
                ./ copiedObj.p (x_J');
            
            % conditioning every single component
            mu    = zeros (length (I), obj.parameters.K);
            sigma = zeros (length (I), length (I), obj.parameters.K);
            prec  = sigma;
            for k = 1:obj.parameters.K
                % per component
                % x_i - mu'_i, with mu'_i = mu_i + sig_IJ * inv(sig_JJ) * (x_J - mu_J)
                mu(:, k) = obj.parameters.mu(I, k) ...
                    + obj.parameters.sigma(I, J, k) ...
                    * inv (obj.parameters.sigma(J, J, k)) ...
                    * (x_J - obj.parameters.mu(J, k));
                sigma(:, :, k) = upperSchurComplement ( ...
                    obj.parameters.sigma(:, :, k), I, J);
                prec(:, :, k) = inv (sigma(:, :, k));
            end % for
            
            obj.parameters.alpha = alpha;
            obj.parameters.mu    = mu;
            obj.parameters.sigma = sigma;
            obj.parameters.prec  = prec;
            
            obj.RVNames = obj.RVNames(~ismember (obj.RVNames, RVNames));
        end % function
        
        % Display model-summary
        function disp (obj)
            if obj.isEstimated
                for k = 1:obj.parameters.K
                    fprintf (1, 'component %d\n', k);
                    disp ('mean');
                    disp (array2table (obj.parameters.mu(:, k)', ...
                        'VariableNames', obj.RVNames));
                    disp ('covariance-matrix');
                    disp (array2table (obj.parameters.sigma(:, :, k), ...
                        'VariableNames', obj.RVNames, ...
                        'RowNames', obj.RVNames));
                    disp ('precision-matrix');
                    disp (array2table (obj.parameters.prec(:, :, k), ...
                        'VariableNames', obj.RVNames, ...
                        'RowNames', obj.RVNames));
                end % for
            else 
                warning ('No parameters has been trained till now.');
            end % if
        end % function
        
        % Returns the probability for given data
        function prob = p (obj, X)
            [M, ~] = size (X);
            
            prob = zeros (M, obj.parameters.K);
            
            for k = 1:obj.parameters.K
                    prob(:, k) = mvnpdf (X, ...
                        obj.parameters.mu(:, k)', obj.parameters.sigma(:, :, k));
            end % for
            prob = sum (prob .* repmat (obj.parameters.alpha', [M, 1]), 2);
        end % function
        
        % Plot the density
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
        
        function inspect (obj)
            inspector (obj);
        end % function
        
        function val = aggregation (obj, type)
            if (nargin < 2)
                type = 'mode'; % maximum of the density
            end % if
            
            if (obj.parameters.K > 1)
                switch (type)
                    case 'mode'
                        relChTol = 1e-6;
                        maxIt    = 50;
                        val      = mean (obj.parameters.mu, 2);
                    
                        it = 1;
                    
                        while (it < maxIt)
                            p = obj.p_perComponent_ (val');
                            
                            s1 = zeros (size (obj.parameters.sigma(:, :, 1)));
                            s2 = zeros (size (obj.parameters.mu(:, 1)));
                            for k = 1:obj.parameters.K
                                tmp = obj.parameters.alpha(k) * p(k) ...
                                    * obj.parameters.prec(:, :, k);
                                s1 = s1 + tmp;
                                s2 = s2 + tmp * obj.parameters.mu(:, k);
                            end % for
                            
                            oldVal = val;
                            val = s1 \ s2;
                            
                            if (abs (val - oldVal) < relChTol)
                                break;
                            end % if
                            
                            it = it + 1;
                        end % while 
                        
                        if (obj.verbose > 0)
                            if (it >= maxIt) 
                                fprintf (1, ...
                                    '[info] Aggregation stoped because the maximum number of iterations has been reached.\n');
                            else
                                fprintf (1, ...
                                    '[info] Aggregation stoped because relative relative-change was below the threshold (it = %d).\n', it);
                            end % if
                        end % if    
                    case 'mean'
                        val = obj.parameters.mu;
                    otherwise
                        error ('No valid aggrigation function: %s.\n', type);
                end % switch
            else
                switch (type)
                    case {'mode', 'mean'}
                        val = obj.parameters.mu;
                    otherwise 
                        error ('No valid aggrigation function: %s.\n', type);
                end % switch
            end % if
        end % function
        
        % Marginalize without side-effects
        function copiedObj = marg_const (obj, RVNames)
            copiedObj = MixtureOfGaussians (obj);
            copiedObj.marg (RVNames);
        end % function
        
        % Conditioning without side-effects
        function copiedObj = cond_const (obj, RVNames, x_J)
            copiedObj = MixtureOfGaussians (obj);
            copiedObj.cond (RVNames, x_J);
        end % function    
        
        % GETTER
        function [alpha, mu, sigma, k] = getDensityParameters (obj)
            alpha = obj.parameters.alpha;
            mu    = obj.parameters.mu;
            sigma = obj.parameters.sigma;
            k     = obj.parameters.K;
        end % function
        
        function RVNames = getRVNames (obj)
            RVNames = obj.RVNames;
        end % function
        
        function mu = getMu (obj)
            mu = obj.parameters.mu;
        end % function
        
        function sigma = getSigma (obj)
            sigma = obj.parameters.sigma;
        end % function
        
    end % methods
    
    methods (Access = 'protected')
        % Estimate model
        function estimateModel_ (obj, data, par)
            [M, ~] = size (data);
            K = par.components;
            if (K < 0)
                error ('Amount of mixture-components must be at least one.');
            end % if
            
            % initialize alphas, mus, sigmas
            obj.parameters.K     = K;
            obj.parameters.alpha = ones (K, 1) * (1 / K);
            obj.parameters.sigma = repmat (cov (data), [1, 1, K]);
            
            if (K > 1)
                obj.parameters.mu    = data(randi (M, 1, K), :)';
                
                % run EM
                obj.EM_ (data, par);
            else
                obj.parameters.mu    = mean (data)';
            end % if
            
            % get precision matrix
            obj.parameters.prec  = zeros (size (obj.parameters.sigma));
            for k = 1:obj.parameters.K
                obj.parameters.prec(:, :, k) = inv (obj.parameters.sigma(:, :, k));
            end % for
        end % function
        
        % EM-Algorithm
        function EM_ (obj, data, par)
            it = 1;
            
            [M, D] = size (data);
            
            alpha = zeros (size (obj.parameters.alpha));
            mu    = zeros (size (obj.parameters.mu));
            sigma = zeros (size (obj.parameters.sigma));
            
            while (it < par.maxIt)
                pData  = obj.p (data);
                
                for k = 1:par.components
                    tmp = obj.parameters.alpha(k) ...
                        * mvnpdf (data, obj.parameters.mu(:, k)', obj.parameters.sigma(:, :, k)) ...
                        ./ pData;
                    
                    % calculate new alpha
                    alpha (k)      = sum (tmp) / M;
                    % calculate new mu
                    mu(:, k)       = sum (repmat (tmp, [1, D]) .* data) / sum (tmp);
                    % calculate new sigma
                    nom = 0;
                    for j = 1:M
                        nom = nom + tmp(j) ...
                            * (data(j, :) - obj.parameters.mu(:, k)')' * (data(j, :) - obj.parameters.mu(:, k)');
                    end % for
                    sigma(:, :, k) = nom / sum (tmp);
                end % for
                
                if (obj.verbose > 1) 
                    fprintf (1, '[deb] sum of alphas = %.8f\n', sum (alpha));
                end % if
                
                % update alpha, mu and sigma
                obj.parameters.alpha = alpha;
                obj.parameters.mu    = mu;
                obj.parameters.sigma = sigma;
                
                % check for relative lh-change
                lhData    = sum (log (pData));
                lhDataNew = sum (log (obj.p (data)));
                if (abs (lhDataNew - lhData) < par.lhTol)
                    break;
                end % if
                
                it = it + 1;
            end % while
            
            if (obj.verbose > 0)
                if (it >= par.maxIt) 
                    fprintf (1, '[info] EM stoped because the maximum number of iterations has been reached.\n');
                else
                    fprintf (1, '[info] EM stoped because relative lh-change was below the threshold (it = %d).\n', it);
                end % if
                fprintf (1, '[info] lh = %.4f\n', lhDataNew);
            end % if            
        end % function
        
        % Get probability for every component
        function prob = p_perComponent_ (obj, X)
            [M, ~] = size (X);
            
            prob = zeros (M, obj.parameters.K);
            
            for k = 1:obj.parameters.K
                    prob(:, k) = mvnpdf (X, ...
                        obj.parameters.mu(:, k)', obj.parameters.sigma(:, :, k));
            end % for
        end % function        
    end % methods 
end % class