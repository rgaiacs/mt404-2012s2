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

% usage: mt404_p03q01();
function mt404_p03q01()
    A = build_matrix();
    % N\'{u}mero de condi\c{c}\~{a}o de A.
    k = cond(A, inf);
    printf('Numero de condicao de $A$: %e\n', k);

    % $b$ como o vetor $e_1$
    b = ones(100, 1);
    x = A \ b;
    % Calculo do res\'{i}duo.
    r = norm(b - A * x, inf);
    printf('Residuo: %e\n', r);
    % Limitante superior para o erro relativo.
    lim_erro = k * r / norm(b, inf);
    printf('Limite superior: %e\n', lim_erro);

    % $b$ tal que a solucao do sistema tenha todas as componentes
    % iguais a 1.
    b = A * ones(100, 1);
    x = A \ b;
    % Calculo do res\'{i}duo.
    r = norm(b - A * x, inf);
    printf('Residuo: %e\n', r);
    % Limitante superior para o erro relativo.
    lim_erro = k * r / norm(b, inf);
    printf('Limite superior: %e\n', lim_erro);
end

% usage: A = build_matrix()
%
% Retorna a matriz a ser utilizada.
function A = build_matrix()
    A = zeros(100);
    for i = 1:99
        A(i,i) = 0.5 + .001 * i;
        A(i,i+1) = -1;
    end
    A(100,100) = .6;
end
