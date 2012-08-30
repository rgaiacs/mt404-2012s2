% Copyright (C) 2012 Raniere Silva <r.gaia.cs@gmail.com>
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% usage: t = mt404_p02();
%
% Benchmark do produto de uma matriz por um vetor.
function t = mt404_p02()
    % valores das dimensoes a serem utilizadas para teste.
    dim = [500, 1000, 2000, 5000, 6000];

    % t(1,:) corresponde ao tempo utilizando [].
    %
    % t(2,:) corresponde ao tempo utilizando sparse.
    %
    % t(3,:) corresponde ao tempo utilizando spdiags.
    t = zeros(3, length(dim));
    for i = 1:length(dim)
        % A \'{e} pentadiagonal
        A = rand(dim(i)) .* sign(conv2(eye(dim(i)),ones(1+2),'same'));
        x = rand(dim(i),1);

        % Sem aproveitamento da estrutura especial.
        tic;
        y = A * x;
        t(1, i) = toc;

        % Utilizando sparse.
        Aux = sparse(A);
        tic;
        y = Aux * x;
        t(2,i) = toc;

        % Utilizando spdiags
        [Aux, d] = spdiags(A);
        tic;
        prod_matrix_spdiags(Aux, x, abs(min(d)), abs(max(d)));
        t(3,i) = toc;
    end
end

% usage: y = prod_matrix_spdiags(A, x, p, q)
%
% Calcula o produto de uma matriz A, cujas diagonais sao armazenadas em
% colunas, por um vetor x.
%
% q eh a amplitudeda banda superior e p eh a amplitude da banda
% inferior.
function y = prod_matrix_spdiags(A, x, p, q)
    % $A$ \'{e} uma matriz banda com banda superior com aplitude $q$ 
    % e banda inferior com amplitude $p$.
    % $x$ \'{e} um vetor
    [m, n] = size(A);
    y = zeros(m,1);
    % i corresponde as linhas
    for i = 1:m
        % j corresponde as colunas
        for j = max(i-p,1):min(i+q,m)
            y(i) = y(i) + spdiags2matrix(A, i, j, p, q) * x(j);
        end
    end
end

% usage: a = spdiags2matrix(B, i, j, p, q)
%
% Retorna o elemento A(i,j) tal que B = spdiags(A).
%
% q eh a amplitudeda banda superior e p eh a amplitude da banda
% inferior da matriz A.
function a = spdiags2matrix(B, i, j, p, q)
    % Get element A(i,j) when B = spdiags(A)
    a = B(j, j - i + p + 1);
end
