% mt404_p01q03b Questao 03b do Projeto 01 de MT404 2012s2
%
% Copyright (C) 2012 Raniere Silva <r.gaia.cs@gmail.com>
% Copyright (C) 2012 Julio Cesar <julioholanda7@gmail.com>
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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
function mt404_p01q03b ()
    n = [2, 5, 10, 20];
    for i = 1:length(n)
        A = triu(rand(n(i)))
        pow2_triu(A)
    end
end

function [ A ] = pow2_triu (A)
    [m, n] = size(A);
    % i corresponde a linhas de A
    for i = 1:m
        % j corresponde a colunas de A
        for j = m:-1:1
            A(i,j) = A(i,:) * A(:,j);
        end
    end
end
