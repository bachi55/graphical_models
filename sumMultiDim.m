%% Function to sum a given matrix over a vector of given dimension
% Input:
%   - A    ... matrix
%   - dims ... dimensions, to sum over
function val = sumMultiDim (A, dims)
    val = A;

    for dim = dims
        val = sum (val, dim);
    end % for
end % function