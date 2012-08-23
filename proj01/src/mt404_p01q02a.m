% mt404_p01q02a Questao 02a do Projeto 01 de MT404 2012s2
%
% Copyright (C) 2012 Raniere Silva <r.gaia.cs@gmail.com>
% Copyright (C) 2012 Julio Cesar <julioholanda7@gmail.com>
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
%
function mt404_p01q02a ()
    n = [2, 5, 10, 20];
    for i = 1:length(n)
        A = rand(n(i));
        x = rand(n(i),1);
        c = meu_prod_mv(A, x);
        gnuo = A * x;
        compara_prod(gnuo, c);
    end
end

function [ c ] = meu_prod_mv (A, x)
    c = zeros(length(x), 1);
    [m, n] = size(A);
    % j corresponde a colunas da matriz A
    for j = 1:n
        % i corresponde a linhas da matriz A
        for i = 1:m
            c(i) = c(i) + A(i,j) * x(j);
        end
    end
end

function compara_prod (m, g)
    ret = m - g;
    ret = abs(ret) >= 10^(-4);
    ret = sum(ret);
    ret = sum(ret);
    if ret == 0
        printf("Matrizes sao iguais.\n");
    else
        printf("Matrizes sao diferentes.\n");
    end
end
