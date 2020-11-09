% Benchmark computer
% benchmark_time = bench;

% Declaration and initialization:
PATH = 'D:\OneDrive - Facultad de Ingeniería\Documents\Academics\fing-aln2020\tarea2';
bcsstk15 = load(fullfile(PATH, 'bcsstk15.mat'));

% Evaluation

% Aca hay que aplicar diferentes estrategias de ordenamiento y acada factorizar una factorizar con chol
% chol(A) Cholesky factorization. https://la.mathworks.com/help/matlab/ref/chol.html

% Analizar la cantidad de memoria utilizada y patrón de la matriz L utilizando las funciones whos y spy.
% spy(S) plots the sparsity pattern of any matrix S. https://la.mathworks.com/help/matlab/ref/spy.html

% Observe la cantidad de elementos distintos de 0 que se obtienen al factorizar una matriz simétrica y definida positiva.

% ----

% ii) Evalúe la resolución de sistemas lineales con el método del gradiente conjugado (pcg en Matlab/Octave) generando sistemas lineales utilizando las matrices de la parte i. 