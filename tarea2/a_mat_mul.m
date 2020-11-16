% Declaration and initialization:
A_2048 = rand(2048);
B_2048 = rand(2048);

A_4096 = rand(4096);
B_4096 = rand(4096);

A_8192 = rand(8192);
B_8192 = rand(8192);

% Evaluation:
disp('Evaluating on 2048')
time_matlab_mul_2048 = timeit(@() matlab_mult(A_2048, B_2048));
time_vector_mul_2048 = timeit(@() vector_mult(A_2048, B_2048));
time_coefs_mul_2048 = timeit(@() coefs_mult(A_2048, B_2048));

disp('    Matlab   Vector    Coefs')
disp([time_matlab_mul_2048, time_vector_mul_2048, time_coefs_mul_2048])

disp('Evaluating on 4096')
time_matlab_mul_4096 = timeit(@() matlab_mult(A_4096, B_4096));
time_vector_mul_4096 = timeit(@() vector_mult(A_4096, B_4096));
time_coefs_mul_4096 = timeit(@() coefs_mult(A_4096, B_4096));

disp('    Matlab   Vector    Coefs')
disp([time_matlab_mul_4096, time_vector_mul_4096, time_coefs_mul_4096])

disp('Evaluating on 8192')
time_matlab_mul_8192 = timeit(@() matlab_mult(A_8192, B_8192));
time_vector_mul_8192 = timeit(@() vector_mult(A_8192, B_8192));
time_coefs_mul_8192 = timeit(@() coefs_mult(A_8192, B_8192));

disp('    Matlab   Vector    Coefs')
disp([time_matlab_mul_8192, time_vector_mul_8192, time_coefs_mul_8192])

% Matrix multiplication methods: (C(m x n) = A(m x p) * B(p x n))

% Using MatLab's tools
function C = matlab_mult(A,B)
    C = A*B;
end

% Multiplying vectors C(i,j) = A(i,:)*B(:,j)
function C = vector_mult(A,B)
    [m,~] = size(A);
    [~,n] = size(B);
    C = zeros(m,n);
    for i = 1:m
       for j = 1:n
           C(i,j) = A(i,:)*B(:,j);
       end
    end
end

% Multiplying coefficients
function C = coefs_mult(A,B)
    [m,p] = size(A);
    [~,n] = size(B);
    C = zeros(m,n);
    for i = 1:m
       for k = 1:p
           for j = 1:n
               C(i,j) = C(i,j) + A(i,k)*B(k,j);
           end
       end
    end
end

