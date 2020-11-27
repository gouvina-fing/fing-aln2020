% A

% Genere 2 matrices completas en Matlab, una de 10x10 y otra de 1024x1024
A_10 = rand(10);
A_1024 = rand(1024);

% Obtenga los valores y vectores propios de dichas matrices utilizando al menos 3 técnicas distintas 
% Evalúe el desempeño de las mismas (tiempo de ejecución, precisión, utilización de memoria, etc.)

% B

% Cargar bcsstk15: https://sparse.tamu.edu/HB/bcsstk15
bcsstk15 = load('bcsstk15.mat');

% Cargar bcsstk01: https://sparse.tamu.edu/HB/bcsstk01
bcsstk01 = load('bcsstk01.mat');

% Repita la parte a) (utilizando técnicas adecuadas para matrices dispersas)

% C

% Leer lena.pgm (greyscale), valores entre 0 y 255
A=imread('lena.pgm');

% Visualizar la imagen
imshow(A)

% C1:
% Utilizando la SVD analice el rango de la imagen

% C2:
% Tome la submatriz correspondiente al primer cuadrante de la imagen
% Repita la parte c1) para dicha submatriz, comparando con el resultado anterior.

% C3:
% Comprima la imagen original usando su SVD, evaluando distintos tamaños de almacenamiento y calidad de imagen.

% C4:
% Proponga una solución para lograr una compresión similar con una mejor calidad de imagen. 