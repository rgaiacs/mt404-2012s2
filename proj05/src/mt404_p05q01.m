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

function mt404_p05q01()
    % Projeto 05, questao 01, de MT404

    % Testes
    test_matrix(1)
    test_matrix(2)
    test_matrix(3)
    test_matrix(4, 10)
    test_matrix(5, 10)
end

function x = solve_t(A, b, is_lower=1)
    % Solve square triangular system of linear equation.
    %
    % :param A: matrix of coeficients.
    %
    % :type A: matrix.
    %
    % :param b: RHS.
    %
    % :type b: vector.
    %
    % :param is_lower: the type of triangular system.
    %
    % :type is_lower: boolean, default 1.
    %
    % :return x: solution for the system.
    %
    % :rtype x: vector.
    dim = issquare(A);
    if dim
        if !is_lower
            % TODO Replace for fliplr and flipud
            A = A(dim:-1:1, dim:-1:1);
            b = b(dim:-1:1);
        end
        % Solve the lower triangular system.
        x = zeros(dim, 1);
        for j = 1:dim
            for i = 1:j - 1
                b(j) = b(j) - A(i, j) * x(j);
            end
            x(j) = b(j) / A(j, j);
        end
        if !is_lower
            x = x(dim:-1:1);
        end
    else
        printf("the input matrix isn't square.");
    end
    
end

function A = build_test_matrix(opt, dim)
    % Build a square matrix to be use as test.
    %
    % :param opt: option to select the matrix to return.
    %
    %       The avaliables values are:
    %
    %           * 1, 2, 3: build pre-defined matrix.
    %           * 4: build 'random' symetric matrix with ``dim`` \times
    %             ``dim``
    %           * 5: build 'random' symetric positive definite matrix with
    %             ``dim \times ``dim``
    % :type opt: integer.
    %
    % :param dim: dimension of the matrix.
    %
    % :type dim: integer.
    %
    % :return: matrix.
    %
    % :rtype: matrix.

    switch(opt)
    case 1
        A = [ 2, -2,   4,  -4;
             -2,  3,  -4,   5;
              4, -4,  10, -10;
             -4,  5, -10,  14];
    case 2
        A = [ 1, -1,  1, -1,  1;
             -1,  2, -2,  2, -2;
              1, -2,  3, -3,  3;
             -1,  2, -3,  4, -4;
              1, -2,  3, -4,  5];
    case 3
        A = [-1, -2,  0,  1;
              2, -3,  2, -1;
              0,  2,  5,  6;
              1, -1,  6, 12];
    case 4
        A = rand(dim);
        A = (A * A) / dim;
    case 5
        A = rand(dim);
        A = (A * A) / dim;
        for i = 1:dim
            if A(i, i) < .1
                A(i, i) = .1;
            end
        end
    otherwise
        A = [];
    end
end

function test_matrix(opt, dim)
    if !exist(dim)
        A = build_test_matrix(opt)
    else
        A = build_test_matrix(opt, dim)
    end

    try
        G = chol(A);
    catch
        printf("can't compute the cholesky factorization of A.");
    end

    if exists(G)
        y = solve_t(G, b, 0);
        x = solve_t(G, y);
    end
end
