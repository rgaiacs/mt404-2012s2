% mt404_p01q01b Questao 01b do Projeto 01 de MT404 2012s2
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
function mt404_p01q01b ()
    n = [3, 3, 4, 4, 6, 6, 10, 10];
    k = [4, 6, 5, 10, 4, 8, 8, 14];
    for i = 1:length(n)
        x = rand(n(i),1);
        y = rand(n(i),1);
        power_prod_ext(x,y,k(i))
    end
end

function [ c ] = power_prod_ext (x, y, k)
    c = (x' * y)^(k - 1) * (x * y');
end
