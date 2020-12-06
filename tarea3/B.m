% Cargar bcsstk15: https://sparse.tamu.edu/HB/bcsstk15
load('bcsstk15.mat');
bcsstk15 = Problem.A;

% Cargar bcsstk01: https://sparse.tamu.edu/HB/bcsstk01
load('bcsstk01.mat');
bcsstk01 = Problem.A;

% https://la.mathworks.com/matlabcentral/answers/104503-how-should-i-compute-the-eigenvectors-of-a-sparse-real-symmetric-matrix

% Repita la parte a) (utilizando técnicas adecuadas para matrices dispersas)

tStart = tic;
[V_k15_eig, D_k15_eig] = eig(full(bcsstk15));
t_k15_eig = toc(tStart);

tStart = tic;
[V_k01_eig, D_k01_eig] = eig(full(bcsstk01));
t_k01_eig = toc(tStart);


tStart = tic;
[V_k15_eig_nb, D_k15_eig_nb] = eig(full(bcsstk15), 'nobalance');
t_k15_eig_nb = toc(tStart);

tStart = tic;
[V_k01_eig_nb, D_k01_eig_nb] = eig(full(bcsstk01), 'nobalance');
t_k01_eig_nb = toc(tStart);


tStart = tic;
[V_k15_eigs, D_k15_eigs] = eigs(bcsstk15, 3948);
t_k15_eigs = toc(tStart);

tStart = tic;
[V_k01_eigs, D_k01_eigs] = eigs(bcsstk01, 48);
t_k01_eigs = toc(tStart);

% Análisis precisión tomando eig como baseline

% Values

error_value_k15_eig_nb = norm(diag(D_k15_eig) - diag(D_k15_eig_nb));
error_value_k15_eigs = norm(diag(D_k15_eig) - diag(D_k15_eigs));

error_value_k01_eig_nb = norm(diag(D_k01_eig) - diag(D_k01_eig_nb));
error_value_k01_eigs = norm(diag(D_k01_eig) - diag(D_k01_eigs));

% Vectors

error_vector_k15_eig_nb = norm(V_k15_eig - V_k15_eig_nb);
error_vector_k15_eigs = norm(V_k15_eig - V_k15_eigs);

error_vector_k01_eig_nb = norm(V_k01_eig - V_k01_eig_nb);
error_vector_k01_eigs = norm(V_k01_eig - V_k01_eigs);