classdef (Abstract) ProbabilityDensityFunction < handle 
    properties (Abstract, Access = 'protected')
        RVNames;
    end % properties
    
    methods (Abstract)    
        % Marginalization
        obj = marg (obj, RVNames);
        % Conditioning
        obj = cond (obj, RVNames, x_J);
        % Display model-summary
        disp (obj);
        % Density values for given data
        p (obj, X);
        % Plot function
        plot (obj, data);
    end % methods
    
    methods (Abstract, Access = 'protected')
        % Estimate model
        estimateModel_ (obj, data, par);
    end % methods 
end % class