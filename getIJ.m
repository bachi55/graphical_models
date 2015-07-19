%% Returns I and J for given RV names. 
% Input:
%   - modelRVNames ... RV names of the model
%   - operRVNames  ... RV names which should be marginalized / contionized
%
% Output:
%   - I ... vector-indicies (according to modelRVNames) of the RV not
%           contained in operRVNames
%   - J ... vector-indicies (according to modelRVNames) of the RV contained 
%           in operRVNames
function [I, J] = getIJ (modelRVNames, operRVNames) 
    [~, J] = ismember (operRVNames, modelRVNames);
    J(J == 0) = [];
    I = 1:length (modelRVNames);
    I(J) = [];
end % function