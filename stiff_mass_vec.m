function [M,S] = stiff_mass_vec(V,T)
% compute the mass and stiff matrix of linear fe on mesh (V,T);
% the local matrix template
nt = size(T,1); 
idx = [1 1 1 2 2 2 3 3 3];
jdx = [1 2 3 1 2 3 1 2 3];

S1 = repmat(0.5*[ 1 -1  0 -1  1  0  0  0  0], nt, 1);
S2 = repmat(0.5*[ 2 -1 -1 -1  0  1 -1  1  0], nt, 1);
S3 = repmat(0.5*[ 1  0 -1  0  0  0 -1  0  1], nt, 1);
S4 = repmat(1/24*[2  1  1 1  2  1 1  1  2], nt, 1);

% comput the local information of all triangles
x21 = V(T(:,2),1) - V(T(:,1),1);
x31 = V(T(:,3),1) - V(T(:,1),1);
y21 = V(T(:,2),2) - V(T(:,1),2);
y31 = V(T(:,3),2) - V(T(:,1),2);
J = x21.*y31 - x31.*y21;           JJ = repmat(J,1,9);
a = (x31.^2 + y31.^2)./J;            aa = repmat(a,1,9);
b = -(x31.*x21 + y31.*y21)./J;   bb = repmat(b,1,9);
c = (x21.^2 + y21.^2)./J;            cc = repmat(c,1,9);

i = T(:, idx);  j = T(:,jdx);
ss = aa.*S1 + bb.*S2 + cc.*S3;
sm = JJ.*S4;

dim = size(V,1);
S = sparse(i, j, ss, dim, dim);
M = sparse(i, j, sm, dim, dim); 