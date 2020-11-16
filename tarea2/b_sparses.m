% Parte I)

% 1. Cargar matriz
load('bcsstk15.mat');
A = Problem.A;

% 2. Variables auxiliares
pct = 100 / numel(A);
nz = nnz(A);
L = chol(A,'lower');
nzl = nnz(L);

% 3. Estrategias de ordenamiento

% 3.1 Column Count - Reordenamiento
p1 = colperm(A);
A1 = A(p1,p1);
nc(1) = nnz(A1);
L1 = chol(A1,'lower');
ncl(1) = nnz(L1);

% 3.2 Reverse Cuthill-McKee - Reordenamiento
p2 = symrcm(A);
A2 = A(p2,p2);
nc(2) = nnz(A2);
L2 = chol(A2,'lower');
ncl(2) = nnz(L2);

% 3.3 Minimum Degree - Reordenamiento
p3 = amd(A);
A3 = A(p3,p3);
nc(3) = nnz(A3);
L3 = chol(A3,'lower');
ncl(3) = nnz(L3);

% 3.4 Nested Dissection Permutation - Reordenamiento
p4 = dissect(A);
A4 = A(p4,p4);
nc(4) = nnz(A4);
L4 = chol(A4,'lower');
ncl(4) = nnz(L4);

% 4. Datos
subplot(2,5,1)
spy(A)
title('A - Matriz Original')
xlabel(sprintf('No nulos = %d (%.2f%%)',nz,nz*pct));

subplot(2,5,2)
spy(A1)
title('A(p1,p1) - Column Count')
xlabel(sprintf('No nulos = %d (%.2f%%)',nc(1),nc(1)*pct));

subplot(2,5,3)
spy(A2)
title('A(p2,p2) - Reverse Cuthill-McKee')
xlabel(sprintf('No nulos = %d (%.2f%%)',nc(2),nc(2)*pct));

subplot(2,5,4)
spy(A3)
title('A(p3,p3) - Minimum Degree')
xlabel(sprintf('No nulos = %d (%.2f%%)',nc(3),nc(3)*pct));

subplot(2,5,5)
spy(A4)
title('A(p4,p4) - Nested Dissection Permutation')
xlabel(sprintf('No nulos = %d (%.2f%%)',nc(4),nc(4)*pct));

subplot(2,5,6)
spy(L)
title('L - Cholesky de A')
xlabel(sprintf('No nulos = %d (%.2f%%)',nzl, nzl*pct));

subplot(2,5,7)
spy(L1)
title('L1 - Cholesky de A1')
xlabel(sprintf('No nulos = %d (%.2f%%)',ncl(1),ncl(1)*pct));

subplot(2,5,8)
spy(L2)
title('L2 - Cholesky de A2')
xlabel(sprintf('No nulos = %d (%.2f%%)',ncl(2),ncl(2)*pct));

subplot(2,5,9)
spy(L3)
title('L3 - Cholesky de A3')
xlabel(sprintf('No nulos = %d (%.2f%%)',ncl(3),ncl(3)*pct));

subplot(2,5,10)
spy(L4)
title('L4 - Cholesky de A4')
xlabel(sprintf('No nulos = %d (%.2f%%)',ncl(4),ncl(4)*pct));

whos

% Parte II)

% 1. Initialization
b_a = rand(size(A,2), 1);         % Real floating-point numbers that are uniformly distributed. Between [0,1]
b_b = randi(10000, size(A,2), 1); % Double integer values drawn from a discrete uniform distribution. Between [0,10000]
b_c = randn(size(A,2), 1);        % Real floating-point numbers that are drawn from a standard normal distribution. With mean 0 and standard deviation 1
tol_p = 1e-14;
tol_np = 1e-14;
maxit = 100000;


% 2. Solve the system of linear equations using Preconditioned Conjugate Gradient
% x      = solution of the system
% flag   = 0 if convergence was succesful.
%        = 1 if pcg iterated maxit iterations but did not converge.
%        = 2 if the preconditioner matrix M or M = M1*M2 is ill conditioned.
%        = 3 if pcg stagnated after two consecutive iterations were the same.
%        = 4 if one of the scalar quantities calculated by the pcg algorithm became too small or too large to continue computing
% relres = relative residual norm(b-A*x)/norm(b). If flag is 0, then relres <= tol.
% iter   = iteration number iter at which x was computed.
% resvec = vector of the residual norm at each iteration, including the first residual norm(b-A*x0).

% 2.0 Using A
[~, flag_p_0_a, relres_p_0_a, iter_p_0_a, ~] = pcg(A, b_a, tol_p, maxit, L, L');
[~, flag_p_0_b, relres_p_0_b, iter_p_0_b, ~] = pcg(A, b_b, tol_p, maxit, L, L');
[~, flag_p_0_c, relres_p_0_c, iter_p_0_c, ~] = pcg(A, b_c, tol_p, maxit, L, L');

[~, flag_np_0_a, relres_np_0_a, iter_np_0_a, ~] = pcg(A, b_a, tol_np, maxit);
[~, flag_np_0_b, relres_np_0_b, iter_np_0_b, ~] = pcg(A, b_b, tol_np, maxit);
[~, flag_np_0_c, relres_np_0_c, iter_np_0_c, ~] = pcg(A, b_c, tol_np, maxit);

% 2.1 Using A: Column Count - Reordenamiento
[~, flag_p_1_a, relres_p_1_a, iter_p_1_a, ~] = pcg(A1, b_a, tol_p, maxit, L1, L1');
[~, flag_p_1_b, relres_p_1_b, iter_p_1_b, ~] = pcg(A1, b_b, tol_p, maxit, L1, L1');
[~, flag_p_1_c, relres_p_1_c, iter_p_1_c, ~] = pcg(A1, b_c, tol_p, maxit, L1, L1');

[~, flag_np_1_a, relres_np_1_a, iter_np_1_a, ~] = pcg(A1, b_a, tol_np, maxit);
[~, flag_np_1_b, relres_np_1_b, iter_np_1_b, ~] = pcg(A1, b_b, tol_np, maxit);
[~, flag_np_1_c, relres_np_1_c, iter_np_1_c, ~] = pcg(A1, b_c, tol_np, maxit);

% 2.2 Using A: Reverse Cuthill-McKee - Reordenamiento
[~, flag_p_2_a, relres_p_2_a, iter_p_2_a, ~] = pcg(A2, b_a, tol_p, maxit, L2, L2');
[~, flag_p_2_b, relres_p_2_b, iter_p_2_b, ~] = pcg(A2, b_b, tol_p, maxit, L2, L2');
[~, flag_p_2_c, relres_p_2_c, iter_p_2_c, ~] = pcg(A2, b_c, tol_p, maxit, L2, L2');

[~, flag_np_2_a, relres_np_2_a, iter_np_2_a, ~] = pcg(A2, b_a, tol_np, maxit);
[~, flag_np_2_b, relres_np_2_b, iter_np_2_b, ~] = pcg(A2, b_b, tol_np, maxit);
[~, flag_np_2_c, relres_np_2_c, iter_np_2_c, ~] = pcg(A2, b_c, tol_np, maxit);

% 2.3 Using A: Minimum Degree - Reordenamiento
[~, flag_p_3_a, relres_p_3_a, iter_p_3_a, ~] = pcg(A3, b_a, tol_p, maxit, L3, L3');
[~, flag_p_3_b, relres_p_3_b, iter_p_3_b, ~] = pcg(A3, b_b, tol_p, maxit, L3, L3');
[~, flag_p_3_c, relres_p_3_c, iter_p_3_c, ~] = pcg(A3, b_c, tol_p, maxit, L3, L3');

[~, flag_np_3_a, relres_np_3_a, iter_np_3_a, ~] = pcg(A3, b_a, tol_np, maxit);
[~, flag_np_3_b, relres_np_3_b, iter_np_3_b, ~] = pcg(A3, b_b, tol_np, maxit);
[~, flag_np_3_c, relres_np_3_c, iter_np_3_c, ~] = pcg(A3, b_c, tol_np, maxit);

% 2.4 Using A: Nested Dissection Permutation - Reordenamiento
[~, flag_p_4_a, relres_p_4_a, iter_p_4_a, ~] = pcg(A4, b_a, tol_p, maxit, L4, L4');
[~, flag_p_4_b, relres_p_4_b, iter_p_4_b, ~] = pcg(A4, b_b, tol_p, maxit, L4, L4');
[~, flag_p_4_c, relres_p_4_c, iter_p_4_c, ~] = pcg(A4, b_c, tol_p, maxit, L4, L4');

[~, flag_np_4_a, relres_np_4_a, iter_np_4_a, ~] = pcg(A4, b_a, tol_np, maxit);
[~, flag_np_4_b, relres_np_4_b, iter_np_4_b, ~] = pcg(A4, b_b, tol_np, maxit);
[~, flag_np_4_c, relres_np_4_c, iter_np_4_c, ~] = pcg(A4, b_c, tol_np, maxit);

% Storing Results
iters_p_a = [iter_p_0_a, iter_p_1_a, iter_p_2_a, iter_p_3_a, iter_p_4_a];
relres_p_a = [relres_p_0_a, relres_p_1_a, relres_p_2_a, relres_p_3_a, relres_p_4_a];
clear iter_p_0_a iter_p_1_a iter_p_2_a iter_p_3_a iter_p_4_a relres_p_0_a relres_p_1_a relres_p_2_a relres_p_3_a relres_p_4_a

iters_p_b = [iter_p_0_b, iter_p_1_b, iter_p_2_b, iter_p_3_b, iter_p_4_b];
relres_p_b = [relres_p_0_b, relres_p_1_b, relres_p_2_b, relres_p_3_b, relres_p_4_b];
clear iter_p_0_b iter_p_1_b iter_p_2_b iter_p_3_b iter_p_4_b relres_p_0_b relres_p_1_b relres_p_2_b relres_p_3_b relres_p_4_b

iters_p_c = [iter_p_0_c, iter_p_1_c, iter_p_2_c, iter_p_3_c, iter_p_4_c];
relres_p_c = [relres_p_0_c, relres_p_1_c, relres_p_2_c, relres_p_3_c, relres_p_4_c];
clear iter_p_0_c iter_p_1_c iter_p_2_c iter_p_3_c iter_p_4_c relres_p_0_c relres_p_1_c relres_p_2_c relres_p_3_c relres_p_4_c

iters_np_a = [iter_np_0_a, iter_np_1_a, iter_np_2_a, iter_np_3_a, iter_np_4_a];
relres_np_a = [relres_np_0_a, relres_np_1_a, relres_np_2_a, relres_np_3_a, relres_np_4_a];
clear iter_np_0_a iter_np_1_a iter_np_2_a iter_np_3_a iter_np_4_a relres_np_0_a relres_np_1_a relres_np_2_a relres_np_3_a relres_np_4_a

iters_np_b = [iter_np_0_b, iter_np_1_b, iter_np_2_b, iter_np_3_b, iter_np_4_b];
relres_np_b = [relres_np_0_b, relres_np_1_b, relres_np_2_b, relres_np_3_b, relres_np_4_b];
clear iter_np_0_b iter_np_1_b iter_np_2_b iter_np_3_b iter_np_4_b relres_np_0_b relres_np_1_b relres_np_2_b relres_np_3_b relres_np_4_b

iters_np_c = [iter_np_0_c, iter_np_1_c, iter_np_2_c, iter_np_3_c, iter_np_4_c];
relres_np_c = [relres_np_0_c, relres_np_1_c, relres_np_2_c, relres_np_3_c, relres_np_4_c];
clear iter_np_0_c iter_np_1_c iter_np_2_c iter_np_3_c iter_np_4_c relres_np_0_c relres_np_1_c relres_np_2_c relres_np_3_c relres_np_4_c