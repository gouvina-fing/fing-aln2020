% Genere 2 matrices completas en Matlab, una de 10x10 y otra de 1024x1024
A_10 = rand(10);
A_1024 = rand(1024);

% Obtenga los valores y vectores propios de dichas matrices utilizando al menos 3 técnicas distintas

% eig:
% http://matlab.izmiran.ru/help/techdoc/ref/eig.html
% [V,D] = eig(A) 
% D: Matriz diagonal de los valores propios
% V matriz cuyas columnas son los vectores propios derechos correspondientes, de manera que A*V = V*D.

% eig is a good, fast, general use eigenvalue/vector solver.
% It is appropriate for use when your matrix is of a realistic size that fits well in memory, and when you need all of the eigenvalues/vectors.
% Sparse matrices do not work at all in eig. Eigs is a solver that is more appropriate for when you need only a limited subset of the eigenvalues/vectors.
% https://xspdf.com/resolution/53574971.html

% Para una matriz real no simétrica A, manteniendo el paso de preliminary balance:
% eig ejecuta las siguientes rutinas de LAPACK:
% DGEEV (with the scaling factor SCLFAC = 2 in DGEBAL, instead of the LAPACK default value of 8)

tStart = tic;
[V_10_eig, D_10_eig] = eig(A_10);
t_10_eig = toc(tStart);

tStart = tic;
[V_1024_eig, D_1024_eig] = eig(A_1024);
t_1024_eig = toc(tStart);

% [V,D] = eig(A,'nobalance') finds eigenvalues and eigenvectors without a preliminary balancing step.
% Ordinarily, balancing improves the conditioning of the input matrix, enabling more accurate computation of the eigenvectors and eigenvalues.
% However, if a matrix contains small elements that are really due to roundoff error, balancing may scale them up to make them as significant as the other elements of the original matrix, leading to incorrect eigenvectors.

% Para una matriz real no simétrica A, sin el paso de preliminary balance:
% eig ejecuta las siguientes rutinas de LAPACK:
% DGEHRD, DORGHR, DHSEQR, DTREVC

tStart = tic;
[V_10_eig_nb, D_10_eig_nb] = eig(A_10, 'nobalance');
t_10_eig_nb = toc(tStart);

tStart = tic;
[V_1024_eig_nb, D_1024_eig_nb] = eig(A_1024, 'nobalance');
t_1024_eig_nb = toc(tStart);


% eigs:
% http://matlab.izmiran.ru/help/techdoc/ref/eigs.html
% Usa un método iterativo encontrar los valores propios. Sirve mejor cuando la matriz es grande, sparse y no simétrica.
% eigs provides the reverse communication required by the Fortran library ARPACK, namely the routines DSAUPD, DSEUPD, DNAUPD, DNEUPD, ZNAUPD, and ZNEUPD.

tStart = tic;
[V_10_eigs, D_10_eigs] = eigs(A_10, 10);
t_10_eigs = toc(tStart);

tStart = tic;
[V_1024_eigs, D_1024_eigs] = eigs(A_1024, 1024);
t_1024_eigs = toc(tStart);

% Análisis precisión tomando eig como baseline

% Values

error_value_10_eig_nb = norm(diag(D_10_eig) - diag(D_10_eig_nb));
error_value_10_eigs = norm(diag(D_10_eig) - diag(D_10_eigs));

error_value_1024_eig_nb = norm(diag(D_1024_eig) - diag(D_1024_eig_nb));
error_value_1024_eigs = norm(diag(D_1024_eig) - diag(D_1024_eigs));

% Vectors

error_vector_10_eig_nb = norm(V_10_eig - V_10_eig_nb);
error_vector_10_eigs = norm(V_10_eig - V_10_eigs);

error_vector_1024_eig_nb = norm(V_1024_eig - V_1024_eig_nb);
error_vector_1024_eigs = norm(V_1024_eig - V_1024_eigs);

% Evalúe el desempeño de las mismas (tiempo de ejecución, precisión, utilización de memoria, etc.)
% con el comando whos se pueden mirar tamaños en memoria.

whos