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
xlabel(sprintf('No nulos = %d',nzl));

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