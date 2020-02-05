function [L, U] = lu_primer(A);
n = size(A,1); 
for k = 1:n-1
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);  
    for i = k+1:n               
        A(i,k+1:end) = A(i,k+1:end) + A(i,k)*A(k,k+1:end);
    end 
end
L = tril(A); U = triu(A) + diag(A);
