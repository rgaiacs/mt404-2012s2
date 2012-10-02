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

function [T, dim, dens, thresh] = mt404_p04(p=0, g=0)
    % Projeto 04 de MT404
    %
    % :param p: 1 to print table.
    %
    % :type p: boolean; default 0.
    %
    % :param g: 1 to generate graphics.
    %
    % :type g: boolean; default 0.
    %
    % :return T: information about the test.
    %
    % :rtype T: list/matrix.
    %
    % :return dim: dimension to be test.
    %
    % :rtype dim: list/vector.
    %
    % :return dens: density to be test.
    %
    % :rtype dens: list/vector.
    %
    % :return thresh: thresh value to be test.
    %
    % :rtype thresh: list/vector.
    %
    % ..  note::
    %
    %     To filter the information about the test by one column, use: ::
    %
    %         > T(T(:, i) == n, :)
    %
    %     To sort the information about the test by one column, use: ::
    %
    %         > [s, s] = sort(T(:, i));
    %         > T(s, :);
    dim = [100, 200, 500, 1000];
    dens = [1, .9, .7, .5, .4, .2, .1];
    thresh = [1, .9, .7, .5, .4, .2, .1, .01, .001];
    % dim = [100, 150, 300];
    % dens = [1, .9, .1];
    % thresh = [.8, .4, .2];
    T = test_lu_thresh(dim, dens, thresh, p, g);
end

function T = test_lu_thresh(dim, dens, thresh, p=0, g=0, v=1)
    % Test LU factorization option THRESH
    %
    % :param dim: dimension to be test.
    %
    % :type dim: list/vector.
    %
    % :param dens: density to be test.
    %
    % :type dens: list/vector.
    %
    % :param thresh: thresh value to be test.
    %
    % :type thresh: list/vector.
    %
    % :param p: 1 to print table.
    %
    % :type p: boolean; default 0.
    %
    % :param g: 1 to generate graphics.
    %
    % :type g: boolean; default 0.
    %
    % :type v: 1 to print the log of execution.
    %
    % :type g: boolean; default 1.
    %
    % :return T: information about the test.
    %
    % :rtype T: list/matrix.
    %
    % An additional input argument THRES, that defines the pivoting
    % threshold can be given.  THRES can be a scalar, in which case it
    % defines UMFPACK pivoting tolerance for both symmetric and
    % unsymmetric cases.  If THRES is a two element vector, then the
    % first element defines the pivoting tolerance for the unsymmetric
    % UMFPACK pivoting strategy and the second the symmetric strategy.
    % By default, the values defined by `spparms' are used and are by
    % default `[0.1, 0.001]'.

    % Disable the figure visualization when create the graphic.
    figure('visible', 'off');

    % Table to save the information.
    % |dimensao    |thresh      |idensidade A |densidade A |densidade L |densidade U |precisao    |tempo       |  
    T = zeros(length(dim) * length(dens) * length(thresh), 8);
    c = 1;
    for n = dim
        for d = dens
            A = sprand(n, n, d);
            s = rand(n, 1);
            b = A * s;
            for t_i = 1:length(thresh)
                if v
                    printf("Testing n = %d, d = %f, thresh = %f\n", n, d, \
                            thresh(t_i));
                    fflush(stdout);
                end

                tic();
                [L, U, P] = lu(A, thresh(t_i));
                y = L \ P * b;
                x = U \ y;
                t = toc();
                erro = norm(x - s, inf) / norm(x, inf);
                
                T(c,:) = [n, thresh(t_i), d, nnz(A) / n^2, nnz(L) / n^2,\
                        nnz(U) / n^2, t, erro];
                c = c + 1;

                % Ilustration.
                if g
                    subplot(length(thresh), 4, t_i * 4 - 3);
                    text(0, .5, sprintf("thresh = %f\ntempo = %f\nerro rel. = %f", thresh(t_i), t, erro))
                    axis([0 1 0 1]);
                    axis("off");
                    subplot(length(thresh), 4, t_i * 4 - 2);
                    axis("on");
                    spy(A);
                    title(sprintf("spy(A), densidade = %f", nnz(A) / n^2));
                    subplot(length(thresh), 4, t_i * 4 - 1);
                    spy(L);
                    title(sprintf("spy(L), densidade = %f", nnz(L) / n^2));
                    subplot(length(thresh), 4, t_i * 4);
                    spy(U);
                    title(sprintf("spy(U), densidade = %f", nnz(U) / n^2));
                end
            end
            if g
                if v
                    printf("Save graphic at %s\n", strcat('../figure/dim', mat2str(n), 'den', mat2str(d, 2), '.png'));
                    fflush(stdout);
                end
                print(strcat('../figure/dim', mat2str(n), 'den', mat2str(d, 2), '.png'), '-dpng' );
                clf;
            end
        end
    end
    if p
        print_table(T);
    end
end

function print_table(T)
    % Print Table.
    %
    % :param T: The table to be print.
    %
    % type T: matrix.
    printf("|dimensao    |thresh      |densidade A |densidade L |densidade U |precisao    |tempo       |\n");  

    printf("|%e|%e|%e|%e|%e|\n", dim, nnz(A)/dim^2,
        nnz(L)/dim^2, nnz(U)/dim^2, thresh(t_i));
end
