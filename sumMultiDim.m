function val = sumMultiDim (A, dims)
    val = A;

    for dim = dims
        val = sum (val, dim);
    end 
end % function