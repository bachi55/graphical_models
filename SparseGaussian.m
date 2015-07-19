classdef SparseGaussian < Gaussian
    methods 
        % Constructor
        function obj = SparseGaussian (dataTable, reg)
            % calls Gaussian constructor
            obj@Gaussian ();
            
            if (nargin < 2)
                reg = 0.5;
            end % if
            
            obj.RVNames = dataTable.Properties.VariableNames;
            data        = table2array (dataTable);
            obj.estimateModel_ (data, struct ('regularization', reg));
            obj.isEstimated = 1;
        end % function
        
        % Plot function
        function plot (obj)
            adjMatrix = (obj.parameters.prec ~= 0);
            view (biograph (adjMatrix - triu (adjMatrix), obj.RVNames, ...
                'ShowArrows', 'off'));
        end % function
    end % methods
    
    methods (Access = 'protected')
        function estimateModel_ (obj, data, par)
            obj.parameters = struct ('mu', mean (data)', 'sigma', cov (data));
            [X, W, ~, ~, iter, ~] = QUIC ('default', obj.parameters.sigma, ...
                par.regularization,  ... % regularization parameter
                1e-6,                ... % convergence threshold
                0,                   ... % verbosity of print statistics
                500);                    % maximum number of Newton iterations to execute
            fprintf (1, '[deb] QUIC converged after %d iterations.\n', iter);
            obj.parameters.sigma = W;
            obj.parameters.prec  = X;
        end % function
    end % methods
end % class