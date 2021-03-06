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

% usage: mt404_p03q02();
function mt404_p03q02()
    for n = [5, 10, 50, 100, 1000]
        H = hilb(n);
        printf('[%d, %d] = size(H)\n', size(H, 1), size(H, 2));

        % N\'{u}mero de condi\c{c}\~{a}o de $H$.
        k = cond(H);
        printf('Numero de condicao de $A$: %e\n', k);

        % $b$ tal que a solucao do sistema tenha todas as componentes
        % iguais a 1.
        b = H * ones(n, 1);
        x = H \ b;
        % Calculo do res\'{i}duo.
        r = norm(b - H * x, inf);
        printf('Residuo: %e\n', r);
        % Limitante superior para o erro relativo.
        lim_erro = k * r / norm(b, inf);
        printf('Limite superior: %e\n', lim_erro);
    end
end
