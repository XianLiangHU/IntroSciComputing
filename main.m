[p,e,t] = initmesh('squareg','Hmax',0.05,'Hgrad',1.99);  %%%%%%%%
V = ([p(1,:)' -p(2,:)']+1)/2;   % original mesh is for domain [-1,1]x[-1,1] 
T = t(1:3,:)'; 

[M,A] = stiff_mass_vec(V,T);

fun_f = @(x,y) 2*pi*pi*sin(pi*x).*sin(pi*y);
fun_u = @(x,y) sin(pi*x).*sin(pi*y);
uu = fun_u(V(:,1),V(:,2));
%trisurf(T,V(:,1),V(:,2),uu);

% Dirichlet boundary on four boundary
bnd = find(abs(V(:,1) < 1e-10 | abs(V(:,1) - 1) < 1e-10 | ...
           abs(V(:,2)) < 1e-10 | abs(V(:,2) - 1) < 1e-10)); 

%% Apply Dirichlet boundary conditions
A(bnd,:) = 0;
A(:,bnd) = 0;
A(bnd,bnd) = speye(length(bnd));
%% the same functional with following for - cycle
% for k = 1:length(bnd_dof)
%     A(bnd_dof(k), bnd_dof(k)) = 1;
% end

% the right hand side
b = M*fun_f(V(:,1),V(:,2));   % <== \int f(x,y)*N_j(x,y) dxdy
b(bnd) = fun_u(V(bnd,1), V(bnd,2));

% and solve it
uh = A \ b;  % bicgstab, cgs

%% by default the first nnode dof is for the nodes
trisurf(T,V(:,1),V(:,2),abs(uh - uu));