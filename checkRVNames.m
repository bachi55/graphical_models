%% Check I and J for given RV names. 
% Input:
%   - modelRVNames ... RV names of the model
%   - operRVNames  ... RV names which should be marginalized / contionized
%
% Output:
%   - isValid ... 1 if all operRVNames are also in modelRVNames, 0
%                 otherwise
function isValid = checkRVNames (modelRVNames, operRVNames) 
    areIn = ismember (operRVNames, modelRVNames);
    
    if nargout > 0
        isValid = all (areIn);
    else 
        if ~all (areIn)
            error ('invalid variable name');
        end % if
    end % if
end % function