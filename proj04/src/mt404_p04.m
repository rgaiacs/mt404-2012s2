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

% usage: mt404_p04();
function mt404_p04()
    % dim = [100, 200, 500, 1000, 2000, 5000];
    % dens = [1, .9, .8, .5, .2, .1, .01, .001];
    % thresh = [1, .9, .8, .5, .4, .2];
    dim = [100];
    dens = [1, .9, .1];
    thresh = [1, .8, .4, .2];
    for n = dim
        test_lu_thresh(n, dens, thresh);
    end
end

function test_lu_thresh(dim, dens, thresh)
    % test_lu_thresh Test LU factorization option THRESH
    %
    % An additional input argument THRES, that defines the pivoting
    % threshold can be given.  THRES can be a scalar, in which case it
    % defines UMFPACK pivoting tolerance for both symmetric and
    % unsymmetric cases.  If THRES is a two element vector, then the
    % first element defines the pivoting tolerance for the unsymmetric
    % UMFPACK pivoting strategy and the second the symmetric strategy.
    % By default, the values defined by `spparms' are used and are by
    % default `[0.1, 0.001]'.

    printf("|dimensao    |densidade A |densidade L |densidade U |thresh      |\n");  
    for d = dens
        A = sprand(dim, dim, d);
        s = rand(dim, 1);
        b = A * s;
        for t = thresh
            [L, U, P] = lu(A, t);
            y = L \ P * b;
            x = U \ y;

            printf("|%e|%e|%e|%e|%e|\n", dim, nnz(A)/dim^2,
                nnz(L)/dim^2, nnz(U)/dim^2, t);
        end
    end
end
