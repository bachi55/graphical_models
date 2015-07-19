%% Function to calculate the upper Schur-Complement
% input: 
%   - a ... matrix A from which the complement schould be calculated
%   - I
%   - J
% output:
%   - S ... upper Schur-Complement of A
function S = upperSchurComplement (A, I, J) 
    S = A(I, I) ...
      - A(I, J) ...
      * inv (A(J, J)) ...
      * A(J, I);
end % function