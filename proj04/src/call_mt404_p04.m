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

mt404_p04('dim', 1)
mt404_p04('zoom', 1, 0, [400], [.5, .25, .1])
% when the matrix size is 400 and the density is 0.01 the matrix is
% singular to machine precision.
mt404_p04('sparse', 1, 0, [400], [.05, .025],
        [1, .75, .5, .25, .1, .05, .025, .01, .005, .0025, .001])
