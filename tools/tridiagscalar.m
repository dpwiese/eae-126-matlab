%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENG 180 Engineering Analysis (Fall 2010)
% Tridiagonal Solver SCALAR
% Daniel Wiese
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves scalar tridiagonal matrix equations. The main diagonals of the matrix are a_diag,
% b_diag, and c_diag. The vector which augments this matrix equation is called d_column. The
% solution to this system is a vector x_sol.
%
% |   b0  c0  0   0   0   ...   0  |   |  x0  |     |  d0  |
% |   a1  b1  c1  0   0   ...   0  |   |  x1  |     |  d1  |
% |   0   a2  b2  c2  0   ...   0  | * |  x2  |  =  |  d2  |
% |   :   :   :   :   :    :    :  |   |  :   |     |  :   |
% |   0   0   0   0 an-2 bn-1 cn-1 |   | xn-1 |     | dn-1 |
% |   0   0   0   0   0   an   bn  |   |  xn  |     |  dn  |
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = tridiagscalar(a, b, c, d)
    a_diag = a;
    b_diag = b;
    c_diag = c;
    d_column = d;

    n_size = length(b_diag);

    % n rows, 1 column
    bb = zeros(n_size,1);
    dd = zeros(n_size,1);

    dd(1) = d_column(1);
    bb(1) = b_diag(1);

    for i = 2:n_size
        bb(i) = b_diag(i)-a_diag(i)*c_diag(i-1)/bb(i-1);
        dd(i) = d_column(i)-a_diag(i)*dd(i-1)/bb(i-1);
    end

    x_sol(n_size) = dd(n_size) / bb(n_size);

    for i = n_size-1:-1:1
        x_sol(i) = (dd(i)-x_sol(i+1)*c_diag(i))/bb(i);
    end

    y = x_sol;
end
