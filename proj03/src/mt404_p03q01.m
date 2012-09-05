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

% usage: mt404_p03q01();
function mt404_p03q01()
    for n = [10, 100]
        A = build_matrix(n);
        printf('[%d, %d] = size(A)\n', size(A, 1), size(A, 2));

        % N\'{u}mero de condi\c{c}\~{a}o de A.
        k = cond(A, inf);
        printf('Numero de condicao de $A$: %e\n', k);

        % $b$ como o vetor de uns
        b = ones(n, 1);
        x = A \ b;
        % Calculo do res\'{i}duo.
        r = norm(b - A * x, inf);
        printf('Residuo: %e\n', r);
        % Limitante superior para o erro relativo.
        lim_erro = k * r / norm(b, inf);
        printf('Limite superior: %e\n', lim_erro);

        % $b$ tal que a solucao do sistema tenha todas as componentes
        % iguais a 1.
        b = A * ones(n, 1);
        x = A \ b;
        % Calculo do res\'{i}duo.
        r = norm(b - A * x, inf);
        printf('Residuo: %e\n', r);
        % Limitante superior para o erro relativo.
        lim_erro = k * r / norm(b, inf);
        printf('Limite superior: %e\n', lim_erro);
    end
end

% usage: A = build_matrix(n)
%
% n eh a dimensao da matriz.
%
% Retorna a matriz a ser utilizada.
function A = build_matrix(n)
    A = zeros(n);
    for i = 1:n-1
        A(i,i) = 0.5 + (.1 / n) * i;
        A(i,i+1) = -1;
    end
    A(n,n) = .6;
end
