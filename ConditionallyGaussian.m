classdef  ConditionallyGaussian < ProbabilityDensityFunction
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
    
    methods
        % Constructor
        % Calls: 
        %   - ConditionallyGaussian():         empty (unestimated) MoG
        %   - ConditionallyGaussian(data):     estimates a MoG from the data
        %   - ConditionallyGaussian(data, par) user defined parameters
        %
        % Input:
        %   - data ... matrix of type Datatable containing the RV names as
        %              column-names
        %   - par  ... struct: 
        %               * 'verbose'    ... output-level: 0=no, 1=info, 2=deb
        function obj = ConditionallyGaussian (varargin)
            % default-constructor
            if (nargin == 0)
                obj.isEstimated = 0;
            % constructor which estimates the density using EM based on the
            % data
            elseif (istable (varargin{1}))
                if (nargin < 2)
                    error ('Please provide the RVNames for the categorical and the continiues variables.');
                end % if
                
                if (nargin < 3)
                    optPar = struct ('verbose', 1); % output-level: 0=no, 1=info, 2=deb
                else
                    optPar = varargin{3};
                end % if
                
                obj.verbose     = optPar.verbose;
                obj.RVNames     = struct ('all'        , {varargin{1}.Properties.VariableNames}, ...
                                          'categorical', {varargin{2}.categorical},              ...
                                          'continuous' , {varargin{2}.continuous}                ...
                                  );
                dataCat         = zeros (size (varargin{1}, 1), numel (obj.RVNames.categorical));
                for i = 1:numel (obj.RVNames.categorical)
                    dataCat(:, i) = categorical ( ...
                                        table2array (varargin{1}(:, i)) ...
                                    );
                end % for
                dataCont        = varargin{1}(:, obj.RVNames.continuous);
                obj.estimateModel_ (dataCat, dataCont, optPar);
                obj.isEstimated = 1;
            % copy-constructor
            elseif (isa (varargin{1}, 'ProbabilityDensityFunction'))
                rhs = varargin{1};
                
                obj.RVNames     = rhs.RVNames;
                obj.indicies    = rhs.indicies;
                obj.infoStr     = rhs.infoStr;
                obj.parameters  = rhs.parameters;
                obj.isEstimated = rhs.isEstimated;
                obj.verbose     = rhs.verbose;
            end % if
        end % ConditionallyGaussian (Constructor)
        
        % Marginalization
        function obj = marg (obj, RVNames)
            if (numel (RVNames) == 0)
                return;
            end % if
            
            checkRVNames (obj.RVNames.all, RVNames);
            
            [ICat,  JCat]  = getIJ (obj.RVNames.categorical, RVNames);
            [ICont, JCont] = getIJ (obj.RVNames.continuous, RVNames);
            
            % calculate weak marginal
            %% p^{-}_J
            p_head_J = squeeze (sumMultiDim (obj.parameters.p, ICat)); 
            
            sizeOfGaussianMatrix = size(obj.parameters.gaussians);
            times = ones (1, ndims(obj.parameters.gaussians));
            times(ICat) = sizeOfGaussianMatrix(ICat);
            p_head_J_origSize = repmat (p_head_J, times);
            
            %% mu^{-}_J
            D = numel (JCont);
            numGaussians = numel (obj.parameters.gaussians);
            
            gaussiansMu = zeros ([numGaussians, D]);
            mu_head_J   = gaussiansMu;
            for i = 1:numGaussians
                gaussiansMu(i, :) = obj.parameters.gaussians(i).marg_const ( ...
                    obj.RVNames.continuous(JCont)).getMu;
                mu_head_J(i, :) = gaussiansMu(i, :) ...
                    * obj.parameters.p(i) / p_head_J_origSize(i);
            end % for
            mu_head_J = reshape (mu_head_J, [sizeOfGaussianMatrix, D]);
            mu_head_J = squeeze (sumMultiDim (mu_head_J, ICat));
            
            mu_head_J_origSize = repmat (mu_head_J, [times, 1]);
            mu_head_J_origSize_lin = reshape (mu_head_J_origSize, [size(gaussiansMu)]);
            
            %% sigma^{-}_JJ
            gaussiansSigma = zeros ([numGaussians, D, D]);
            sigma_head_JJ  = gaussiansSigma;
            for i = 1:numGaussians
                gaussiansSigma(i, :, :) = obj.parameters.gaussians(i).marg_const ( ...
                    obj.RVNames.continuous(JCont)).getSigma;
                tmp = gaussiansMu(i, :) - mu_head_J_origSize_lin(i, :);
                sigma_head_JJ(i, :, :) = (tmp' * tmp + squeeze (gaussiansSigma(i, :, :))) ...
                                       * obj.parameters.p(i) / p_head_J_origSize(i);
            end % for
            sigma_head_JJ = reshape (sigma_head_JJ, [sizeOfGaussianMatrix, D, D]);
            sigma_head_JJ = squeeze (sumMultiDim (sigma_head_JJ, ICat));
            
            %% update the obj
            obj.parameters.p         = p_head_J;
            obj.parameters.gaussians = repmat (Gaussian2, size (p_head_J));
            
            
                  
            obj.RVNames.all        = RVNames;
            obj.RVNames.categorical = obj.RVNames.categorical(JCat);
            obj.RVNames.continuous = obj.RVNames.continuous(JCont);
        end % function
        % Conditioning
        function obj = cond (obj, RVNames, x_J)
        end % function
        % Display model-summary
        function disp (obj)
        end
        % Density values for given data
        function p (obj, X)
        end 
        % Plot function
        function plot (obj, data)
        end 
        
        % GETTER
        function [p, gaussians] = getDensityParameters (obj)
            p         = obj.parameters.p;
            gaussians = obj.parameters.gaussians;
        end % getDensityParameters
    end % methods (public)
    
    methods (Access = 'protected')
        function estimateModel_ (obj, dataCat, dataCont, optPar)
            numcategoricalRV = size (dataCat,  2);
            numContonuousRV = size (dataCont, 2);
            
            if (numcategoricalRV > 1)
                pSize           = zeros (1, numcategoricalRV);
                for i = 1:numcategoricalRV
                    pSize(i) = numel (unique (dataCat(:, i)));
                end % for 
                obj.parameters.p         = zeros (pSize);
                obj.parameters.gaussians = repmat (Gaussian2, pSize);
                
                obj.parameters.mus       = zeros ([pSize, numContonuousRV]);
                obj.parameters.sigma     = zeros ( ...
                    [pSize, numContonuousRV, numContonuousRV]);
            else 
                numCategories = numel (unique (dataCat));
                obj.parameters.p         = zeros (numCategories, 1);
                obj.parameters.gaussians = repmat (Gaussian2, [numCategories, 1]);
                
                obj.parameters.mus       = zeros ([numCategories, numContonuousRV]);
                obj.parameters.sigma     = zeros ( ...
                    [numCategories, numContonuousRV, numContonuousRV]);
            end % if
            
            unCatDataRows = unique (dataCat, 'rows');
            
            % maximize L(p, mu, sigma) = L(p) * PROD L(mu, sigma)
            for i = 1:size (unCatDataRows, 1)
                unCatDataRow = unCatDataRows(i, :);
                
                % x in Omega_x with x = unCatDataRow
                ismem = ismember (dataCat, unCatDataRow, 'rows');
                
                % maximize L(p)
                obj.parameters.p(double (unCatDataRow)) = ...
                    sum (ismem) / size (dataCat, 1);
                
                % maximize L(mu, sigma)
                obj.parameters.gaussians(double (unCatDataRow)) = ...
                    Gaussian2 (dataCont(ismem, :), optPar);
                
                str = strcat (sprintf('%d,', double (unCatDataRow)), ':');
                eval ( ...
                    ['obj.parameters.mus(' str ') = transp (mean (dataCont{ismem, :}))']);
                
                str = strcat (sprintf('%d,', double (unCatDataRow)), ':,:');
                eval ( ...
                    ['obj.parameters.sigma(' str ') = cov (dataCont{ismem, :})']);
            end % for            
        end % estimateModel_ 
    end % methods (protected)
end % classdef