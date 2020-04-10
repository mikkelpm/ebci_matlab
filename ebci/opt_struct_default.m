function opt_struct = opt_struct_default()

    % Struct with default optimization options

    opt_struct = struct;
    opt_struct.fzero = optimset('Display', 'off'); % Default options for "fzero"
    opt_struct.fminbnd = optimset('Display', 'off'); % Default options for "fminbnd"
    opt_struct.linprog = optimoptions('linprog', 'Display', 'off'); % Default options for "linprog"
    opt_struct.check = true; % Double-check critical value calculation?
    opt_struct.numgrid = 5e3; % Number of grid points for linear program that computes maximal non-coverage

end