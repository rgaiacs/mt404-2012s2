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

function mt404_p04(f_name='dim', to_print=0, to_group=1, dim=[100, 200, 400], dens=[1, .75, .5, .25, .1], thresh=[1, .75, .5, .25, .1])
    % Projeto 04 de MT404
    %
    % :param f_name: name of the file to save the information.
    %
    % :type f_name: string; default 'benchmark.csv'.
    %
    % :param to_print: 1 to generate graphics.
    %
    % :type to_print: boolean; default 0.
    %
    % :param to_group: 1 to group graphics of same dimension.
    %
    % :type to_group: boolean; default 1.
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
    % ..  note::
    %
    %     The information is save in a cvs file.
    T = test_lu_thresh(dim, dens, thresh, to_print, to_group, f_name);
    printf("Save information to %s\n", strcat(f_name, '.csv'));
    fflush(stdout);
    csvwrite(strcat(f_name, '.csv'), T);
end

function T = test_lu_thresh(dim, dens, thresh, to_print=0, to_group=1, f_name='dim', v=1)
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
    % :param to_print: 1 to generate graphics.
    %
    % :type to_print: boolean; default 0.
    %
    % :param to_group: 1 to group graphics of same dimension.
    %
    % :type to_group: boolean; default 1.
    %
    % :param f_name: prefix of the name of figures to save.
    %
    % :type f_name: string; default ''.
    %
    % :param v: 1 to print the log of execution.
    %
    % :type v: boolean; default 1.
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
        for d_i = 1:length(dens)
            A = sprand(n, n, dens(d_i));
            s = rand(n, 1);
            b = A * s;
            for t_i = 1:length(thresh)
                if v
                    printf("Testing n = %d, d = %f, thresh = %f\n", n, dens(d_i), \
                            thresh(t_i));
                    fflush(stdout);
                end

                tic();
                [L, U, P] = lu(A, thresh(t_i));
                y = L \ P * b;
                x = U \ y;
                t = toc();
                erro = norm(x - s, inf) / norm(x, inf);
                
                T(c,:) = [n, thresh(t_i), dens(d_i), nnz(A) / n^2, nnz(L) / n^2,\
                        nnz(U) / n^2, t, erro];
                c = c + 1;

                % Ilustration.
                if to_print
                    if to_group
                        subplot(length(dens), length(thresh), d_i * length(dens) - length(thresh) + t_i);
                        hold on;
                        spy(L, 'sb');
                        spy(U, 'sb');
                        spy(A, 'sm');
                        title(sprintf("densidade = %f\nthresh = %f", nnz(A) / n^2, thresh(t_i)));
                    else
                        if v
                            printf("Save graphic at %s\n",
                            strcat(f_name, mat2str(n), '-', mat2str(d_i), '-', mat2str(t_i), '.png'));
                            fflush(stdout);
                        end
                        hold on;
                        spy(L, 'sb');
                        spy(U, 'sb');
                        spy(A, 'sm');
                        title(sprintf("densidade = %f\nthresh = %f", nnz(A) / n^2, thresh(t_i)));
                        print(strcat(f_name, mat2str(n), '-',
                        mat2str(d_i), '-', mat2str(t_i), '.png'), '-dpng', '-S800,800');
                        clf;
                    end
                end
            end
        end
        if to_print && to_group
            if v
                printf("Save graphic at %s\n", strcat(f_name, mat2str(n), '.png'));
                fflush(stdout);
            end
            print(strcat(f_name, mat2str(n), '.png'), '-dpng', '-S3000,3000');
            clf;
        end
    end
end
