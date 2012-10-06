% Copyright (C) 2012 Raniere Silva <r.gaia.cs@gmail.com>
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function mt404_p05q02()
    % Projeto 05, questao 02, de MT404
    n = 0;
    is_singular = 0;
    while !is_singular
        n = n + 1;

        % -- Function File:  hilb (N)
        %     Return the Hilbert matrix of order N.  The i, j element of a
        %     Hilbert matrix is defined as
        %
        %          H (i, j) = 1 / (i + j - 1)
        H = hilb(n);
        try
            chol(H);
        catch
            is_singular = 1;
        end
    end

    printf('The high order Hilbert matrix that is symmetric positive definite is %d\n', n);
end
