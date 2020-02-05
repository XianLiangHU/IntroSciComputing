function x = gauss(A, b);
n = size(A,1);
for k = 1:n-1
    A(k,:) = A(k,:)/A(k,k);          
    for j = k+1:n              
        factor = -A(j,k)/A(k,k);   
        A(j,k:end) = A(j,k:end) + factor*A(k,k:end);
        b(j) = b(j) + factor*b(k);
    end 
end
b(n) = b(n)/A(n,n);      
for k = n:-1:2
	b(1:k-1) = b(1:k-1) - A(1:k-1,k)*b(k); 
end
x = b;         
