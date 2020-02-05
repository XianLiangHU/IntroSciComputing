load clown.mat; 
[U,S,V]=svd(X); 
image(U(:,1:k)*S(1:k,1:k)*V(:,1:k)');
