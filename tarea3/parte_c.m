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

lena_192 = lena_U(:,1:192)*lena_S(1:192,1:192)*lena_V(:,1:192)';
lena_128 = lena_U(:,1:128)*lena_S(1:128,1:128)*lena_V(:,1:128)';
lena_96 = lena_U(:,1:96)*lena_S(1:96,1:96)*lena_V(:,1:96)';
lena_64 = lena_U(:,1:64)*lena_S(1:64,1:64)*lena_V(:,1:64)';
lena_32 = lena_U(:,1:32)*lena_S(1:32,1:32)*lena_V(:,1:32)';
lena_16 = lena_U(:,1:16)*lena_S(1:16,1:16)*lena_V(:,1:16)';
lena_8 = lena_U(:,1:8)*lena_S(1:8,1:8)*lena_V(:,1:8)';

imwrite(lena_192, 'tarea3/results/lena_192.pgm')
imwrite(lena_128, 'tarea3/results/lena_128.pgm')
imwrite(lena_96, 'tarea3/results/lena_96.pgm')
imwrite(lena_64, 'tarea3/results/lena_64.pgm')
imwrite(lena_32, 'tarea3/results/lena_32.pgm')
imwrite(lena_16, 'tarea3/results/lena_16.pgm')
imwrite(lena_8, 'tarea3/results/lena_8.pgm')

[lena_psnr_192, lena_snr_192] = psnr(lena, lena_192);
lena_rmse_192 = sqrt(immse(lena, lena_192));
[lena_psnr_128, lena_snr_128] = psnr(lena, lena_128);
lena_rmse_128 = sqrt(immse(lena, lena_128));
[lena_psnr_96, lena_snr_96] = psnr(lena, lena_96);
lena_rmse_96 = sqrt(immse(lena, lena_96));
[lena_psnr_64, lena_snr_64] = psnr(lena, lena_64);
lena_rmse_64 = sqrt(immse(lena, lena_64));
[lena_psnr_32, lena_snr_32] = psnr(lena, lena_32);
lena_rmse_32 = sqrt(immse(lena, lena_32));
[lena_psnr_16, lena_snr_16] = psnr(lena, lena_16);
lena_rmse_16 = sqrt(immse(lena, lena_16));
[lena_psnr_8, lena_snr_8] = psnr(lena, lena_8);
lena_rmse_8 = sqrt(immse(lena, lena_8));

lena_size = 256^2;
lena_svd_size_192 = 2*(192^2) + 192;
lena_svd_size_128 = 2*(128^2) + 128;
lena_svd_size_96 = 2*(96^2) + 96;
lena_svd_size_64 = 2*(64^2) + 64;
lena_svd_size_32 = 2*(32^2) + 32;
lena_svd_size_16 = 2*(16^2) + 16;
lena_svd_size_8 = 2*(8^2) + 8;

whos

% SVD Compression, image size: lena_size elements

figure;

subplot(2,3,1)
imshow(lena_192)
title('Compressed, k = 192')
xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.4f dB; SNR: %.4f db; RMSE: %.4f', lena_svd_size_192, 100*lena_size/lena_svd_size_192, lena_psnr_192, lena_snr_192, lena_rmse_192));

subplot(2,3,2)
imshow(lena_128)
title('Compressed, k = 128')
xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.4f dB; SNR: %.4f dB; RMSE: %.4f', lena_svd_size_128, 100*lena_size/lena_svd_size_128, lena_psnr_128, lena_snr_128, lena_rmse_128));

subplot(2,3,3)
imshow(lena_96)
title('Compressed, k = 96')
xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.4f dB; SNR: %.4f dB; RMSE: %.4f', lena_svd_size_96, 100*lena_size/lena_svd_size_96, lena_psnr_96, lena_snr_96, lena_rmse_96));

subplot(2,3,4)
imshow(lena_64)
title('Compressed, k = 64')
xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.4f dB; SNR: %.4f dB; RMSE: %.4f', lena_svd_size_64, 100*lena_size/lena_svd_size_64, lena_psnr_64, lena_snr_64, lena_rmse_64));

subplot(2,3,5)
imshow(lena_32)
title('Compressed, k = 32')
xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.4f dB; SNR: %.4f dB; RMSE: %.4f', lena_svd_size_32, 100*lena_size/lena_svd_size_32, lena_psnr_32, lena_snr_32, lena_rmse_32));

subplot(2,3,6)
imshow(lena_16)
title('Compressed, k = 16')
xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.4f dB; SNR: %.4f dB; RMSE: %.4f', lena_svd_size_16, 100*lena_size/lena_svd_size_16, lena_psnr_16, lena_snr_16, lena_rmse_16));

% subplot(2,3,7)
% imshow(lena_8)
% title('Compressed, k = 8')
% xlabel(sprintf('%d elements svd (%.2f%% compression ratio) \n PSNR: %.2f dB; SNR: %.2f dB; MSE: %.2f', lena_svd_size_8, 100*lena_size/lena_svd_size_8, lena_psnr_8, lena_snr_8, lena_rmse_8));

% C4:
% Proponga una solución para lograr una compresión similar con una mejor calidad de imagen.
