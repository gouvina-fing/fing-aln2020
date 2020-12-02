% A

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

[V_10_eig, D_10_eig] = eig(A_10);
[V_1024_eig, D_1024_eig] = eig(A_1024);

% [V,D] = eig(A,'nobalance') finds eigenvalues and eigenvectors without a preliminary balancing step.
% Ordinarily, balancing improves the conditioning of the input matrix, enabling more accurate computation of the eigenvectors and eigenvalues.
% However, if a matrix contains small elements that are really due to roundoff error, balancing may scale them up to make them as significant as the other elements of the original matrix, leading to incorrect eigenvectors.

% Para una matriz real no simétrica A, sin el paso de preliminary balance:
% eig ejecuta las siguientes rutinas de LAPACK:
% DGEHRD, DORGHR, DHSEQR, DTREVC

[V_10_eig_nb, D_10_eig_nb] = eig(A_10, 'nobalance');
[V_1024_eig_nb, D_1024_eig_nb] = eig(A_1024, 'nobalance');

% eigs:
% http://matlab.izmiran.ru/help/techdoc/ref/eigs.html
% Usa un método iterativo encontrar los valores propios. Sirve mejor cuando la matriz es grande, sparse y no simétrica.
% eigs provides the reverse communication required by the Fortran library ARPACK, namely the routines DSAUPD, DSEUPD, DNAUPD, DNEUPD, ZNAUPD, and ZNEUPD.

[V_10_eigs, D_10_eigs] = eigs(A_10, 10);
[V_1024_eigs, D_1024_eigs] = eigs(A_1024, 1024);

% Estos metodos no andan? Borrar?
% l = ww(A_10, 0.00001);
% [e_val,e_vec] = Power_method_dominant_eig(A_10)
% [Q, R] = QR_GramSchmidt(A_10)

% Evalúe el desempeño de las mismas (tiempo de ejecución, precisión, utilización de memoria, etc.)
% con el comando whos se pueden mirar tamaños en memoria.

% B

% Cargar bcsstk15: https://sparse.tamu.edu/HB/bcsstk15
load('bcsstk15.mat');
bcsstk15 = Problem.A;

% Cargar bcsstk01: https://sparse.tamu.edu/HB/bcsstk01
load('bcsstk01.mat');
bcsstk01 = Problem.A;

% https://la.mathworks.com/matlabcentral/answers/104503-how-should-i-compute-the-eigenvectors-of-a-sparse-real-symmetric-matrix

% Repita la parte a) (utilizando técnicas adecuadas para matrices dispersas)

% tic; 
[V_k15_eigs, D_k15_eigs] = eigs(bcsstk15, 3948);
[V_k15_eig, D_k15_eig] = eig(full(bcsstk15));

[V_k01_eigs, D_k01_eigs] = eigs(bcsstk01, 48);
[V_k01_eig, D_k01_eig] = eig(full(bcsstk01));

% C

% Leer lena.pgm (greyscale), valores entre 0 y 255
lena = imread('lena.pgm');
lena = im2double(lena);

% Visualizar la imagen
imshow(lena)

% C1:
% Utilizando la SVD analice el rango de la imagen
lena_svd = svd(lena);
[lena_U, lena_S, lena_V] = svd(lena); % performs a singular value decomposition of matrix A, such that A = U*S*V'.

% rang(A) = cantidad de valores singulares de A distintos de 0.
lena_rango = nnz(lena_svd);

% C2:
% Tome la submatriz correspondiente al primer cuadrante de la imagen
[rows, columns] = size(lena);
lena_1st_quarter = lena(1:rows/2, 1:columns/2);
imshow(lena_1st_quarter)

% Repita la parte c1) para dicha submatriz, comparando con el resultado anterior.
lena_1q_svd = svd(lena_1st_quarter);
[lena_1q_U, lena_1q_S, lena_1q_V] = svd(lena_1st_quarter); % performs a singular value decomposition of matrix A, such that A = U*S*V'.

% rang(A) = cantidad de valores singulares de A distintos de 0.
lena_1q_rango = nnz(lena_1q_svd);

% C3:
% Comprima la imagen original usando su SVD, evaluando distintos tamaños de almacenamiento y calidad de imagen.

% Descomponer la matriz de imagen y luego comprimir la información utilizando solo algunos valores singulares dependiendo de la calidad que deseamos obtener.
% Para recuperar la matriz original después de aplicarle la SVD podemos utilizar la definición de suma de tripletas.
% Para comprimir podemos reconsturirla hasta un numero de valore singulares
% entre 1 y el rango de la matriz (1 < k < r)

lena_256 = lena_U*lena_S*lena_V';
lena_128 = lena_U(:,1:128)*lena_S(1:128,1:128)*lena_V(:,1:128)';
lena_64 = lena_U(:,1:64)*lena_S(1:64,1:64)*lena_V(:,1:64)';
lena_32 = lena_U(:,1:32)*lena_S(1:32,1:32)*lena_V(:,1:32)';
lena_16 = lena_U(:,1:16)*lena_S(1:16,1:16)*lena_V(:,1:16)';
lena_8 = lena_U(:,1:8)*lena_S(1:8,1:8)*lena_V(:,1:8)';

whos

figure;

subplot(2,3,1)
imshow(lena_256)
title('Original, k = 256')

subplot(2,3,2)
imshow(lena_128)
title('Compressed, k = 128')

subplot(2,3,3)
imshow(lena_64)
title('Compressed, k = 64')

subplot(2,3,4)
imshow(lena_32)
title('Compressed, k = 32')

subplot(2,3,5)
imshow(lena_16)
title('Compressed, k = 16')

subplot(2,3,6)
imshow(lena_8)
title('Compressed, k = 8')

% C4:
% Proponga una solución para lograr una compresión similar con una mejor calidad de imagen. 

function [Q R]=QR_GramSchmidt(A)
[n n]=size(A);
Q=zeros(n);
R=zeros(n);
R(1,1)=norm(A(:,1));
Q(:,1)=A(:,1)/R(1,1);

for j=2:n
    Q(:,j)=A(:,j);
    for i=1:j-1
    Q(:,j)=Q(:,j)-A(:,j)'*Q(:,i)*Q(:,i);

    R(i,j)=A(:,j)'*Q(:,i);

    fprintf(['At iteration i = ',num2str(i),' QR is equal to ']) % print Ai
    Q*R
end
R(j,j)=norm(Q(:,j));
Q(:,j)=Q(:,j)/norm(Q(:,j));
fprintf(['At iteration j = ',num2str(j),' QR is equal to '])
    Q*R
end   
end

function l = ww(A,E)
    x = A(:,1);
    l = 0;
    blad = E;        % starting value of error
    while blad>=E
        y = A*x;
        blad = l;
        l = x.*y;  % Rayleigh
        m = x.*x;
        l = l/m;           % eigenvalue
        blad = abs(l - blad);         % error
        x = y;
    end
end

function [e_val,e_vec] = Power_method_dominant_eig(A,x0,error,no_itr)
    if nargin < 3 , error = 0.00001 ;no_itr = 20;  end
    if nargin < 2 , error = 0.00001 ;no_itr = 20; x0=ones(size(A));x0=x0(:,1);  end
    if nargin < 4 , error = 0.00001 ; end
    x1 = A*x0;
    x1 = x1/max(x1);
    norm1 =norm(x1-x0);
    i = 0;
    while norm1 > 0.01
        x2 = A*x1;
        x2 = x2/max(real(x2));
        norm1 = norm(x1-x2);
        x1 = x2;
        i = i + 1;
        if no_itr > i && norm1 < error , break , end
    end
    e_vec = x1;
    e_val = ((A*x1)'*x1)/(x1'*x1);
    
end