classdef Gaussian2 < MixtureOfGaussians
    methods 
        % Constructor
        % Calls: 
        %   - Gaussian2():         empty (unestimated) MoG
        %   - Gaussian2(data):     estimates a MoG from the data
        %   - Gaussian2(data, par) user defined parameters
        %
        % Input:
        %   - data ... matrix of type Datatable containing the RV names as
        %              column-names
        %   - par  ... struct: 
        %               * 'verbose'    ... output-level: 0=no, 1=info, 2=deb
        function obj = Gaussian2 (varargin)
            if (nargin == 2)
                varargin{2}.components = 1;
            end 
            
            obj@MixtureOfGaussians (varargin{:});
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
        
        % Marginalize without side-effects
        function copiedObj = marg_const_ (obj, RVNames)
            copiedObj = Gaussian2 (obj);
            copiedObj.marg (RVNames);
        end % function
        
        function scaleMu (obj, factor)
            obj.parameters.mu = factor * obj.parameters.mu;
        end % function
        
        % GETTER
        function [mu, sigma] = getDensityParameters (obj)
            mu    = obj.parameters.mu;
            sigma = obj.parameters.sigma;
        end % function
    end % methods
end % classdef