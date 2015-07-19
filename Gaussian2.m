classdef Gaussian2 < MixtureOfGaussians
    methods 
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
        
        % SETTER
        
        % OPERATOR
%         function outObj = plus (obj1, obj2)       
%         end % function
    end % methods
end % classdef