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

function mt404_p05q01(show=0, f_name='mt404_p05q01.csv')
    % Projeto 05, questao 01, de MT404

    % Testes
    info = [];
    [erro, res, opt, dim] = test_matrix(1);
    info = [info; [1, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(2);
    info = [info; [2, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(3);
    info = [info; [3, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(4, 10);
    info = [info; [4, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(4, 50);
    info = [info; [5, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(4, 100);
    info = [info; [6, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(5, 10);
    info = [info; [7, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(5, 50);
    info = [info; [8, opt, dim, erro, res]];
    [erro, res, opt, dim] = test_matrix(5, 100);
    info = [info; [9, opt, dim, erro, res]];
    if show
        disp(info)
    end
    printf(cstrcat('writing information at ', f_name, '.\n'));
    csvwrite(f_name, info);
end

function x = solve_with_chol(A, b)
    % Try to solve a system of linear equation, where $A$ is symmetric, using
    % the Cholesky Decomposition.
    %
    % :param A: matrix of coeficients.
    %
    % :type A: matrix.
    %
    % :param b: RHS.
    %
    % :type b: vector.
    %
    % :return x: solution for the system.
    %
    % :rtype x: vector.
    try
        % `chol' is a function from the file
        % /usr/lib/octave/3.2.4/oct/i686-pc-linux-gnu/chol.oct 
        %
        % -- Loadable Function: R = chol (A)
        % -- Loadable Function: [L, ...] = chol (..., 'lower')
        % Compute the Cholesky factor, R, of the symmetric positive
        % definite matrix A, where
        %
        %     R' * R = A.
        %
        % Called with either a sparse or full matrix and using the 'lower'
        % flag, `chol' returns the lower triangular factorization such that
        %
        %     L * L' = A.
        %
        % In general the lower triangular factorization is significantly faster
        % for sparse matrices.
        L = chol(A, 'lower');
    catch
        printf("can't compute the cholesky factorization of A whos dimension are %d \\times %d.\n", size(A, 1), size(A, 2));
        f_name = strcat(mat2str(floor(1000 * rand())), '.dat');
    end

    if exist('L', 'var')
        y = solve_tl(L, b);
        x = solve_tu(L', y);
    else
        x = [];
    end
end

function x = solve_tl(A, b)
    % Solve square lower triangular system of linear equation.
    %
    % :param A: lower triangular matrix of coeficients.
    %
    % :type A: matrix.
    %
    % :param b: RHS.
    %
    % :type b: vector.
    %
    % :return x: solution for the system.
    %
    % :rtype x: vector.
    dim = issquare(A);
    if dim
        x = zeros(dim, 1);
        for j = 1:dim
            for i = 1:j - 1
                b(j) = b(j) - A(i, j) * x(j);
            end
            x(j) = b(j) / A(j, j);
        end
    else
        printf("the input matrix isn't square.");
    end
end

function x = solve_tu(A, b)
    % Solve square upper triangular system of linear equation.
    %
    % :param A: matrix of coeficients.
    %
    % :type A: matrix.
    %
    % :param b: RHS.
    %
    % :type b: vector.
    %
    % :return x: solution for the system.
    %
    % :rtype x: vector.
    dim = issquare(A);
    if dim
        x = zeros(dim, 1);
        for j = dim:-1:1
            for i = dim:-1:j + 1
                b(j) = b(j) - A(i, j) * x(j);
            end
            x(j) = b(j) / A(j, j);
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
        A = tril(rand(dim));
        A = (A * A');
    case 5
        A = eye(dim) + tril(rand(dim));
        A = (A * A');
    otherwise
        A = [];
    end
end

function [erro, res, opt, dim] = test_matrix(opt, dim)
    if !exist('dim', 'var')
        A = build_test_matrix(opt);
        dim = size(A, 1);
    else
        A = build_test_matrix(opt, dim);
    end
    s = rand(size(A, 2), 1);
    b = A * s;

    x = solve_with_chol(A, b);
    if x
        erro = norm(x - s, inf) / norm(x, inf);
        res = norm(b - A * x, inf) / norm(b, inf);
    else
        erro = inf;
        res = inf;
    end
end
