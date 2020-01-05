% This is a demonstration for solving elliptic bvp with mFEM
% 

%% 0. the finite element mesh  %%%%%%%%
t = (0:0.05:1)'; 
[V, T] = msh_init_rect(t, t, 'type4B');  

% another mesh: delaunay type
% [p,e,t] = initmesh('squareg','Hmax',0.15,'Hgrad',1.99);  %%%%%%%%
% V = ([p(1,:)' -p(2,:)']+1)/2;   % original mesh is for domain [-1,1]x[-1,1] 
% T = t(1:3,:)'; 

%   build the finite element mesh
[Nodes, Edges, Elements] = msh_build_fem(V, T);

% 1.2 specify  the Dirichlet boundary on four boundary
edge_center = (Nodes(Edges(:,1),:) + Nodes(Edges(:,2),:))/2;
DirichletSides = find(abs(edge_center(:,1)) < 1e-10 | abs(edge_center(:,1) - 1) < 1e-10 | ...
                      abs(edge_center(:,2)) < 1e-10 | abs(edge_center(:,2) - 1) < 1e-10); 

% 1.3 specify the nuemann boundary
NeumannSides = [];  % 2 x nn matrix, each column means a edge

%% 2. the finite element information
%  eval the basis functions' value at the quadrature points
n_quad = 7;     %%%%%%%%%
[quad_weight, quad_pnt] = quad_rule(n_quad);
FE_Order = 2;   %%%%%%%%%
[bas_val, bas_grad_xi, bas_grad_eta] = fe_lagrange(quad_pnt(:,1), quad_pnt(:,2), FE_Order);
[dof_map_elem, dof_map_edge, dof_xyz, dof_patten_side, dof_xieta_elem, n_dof] = dof_lagrange(Nodes, Edges, Elements, FE_Order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. element-wise stiff matrix: according to different type of PDEs
% % 1.1 parameter functions in the elliptic equation
fun_D = @(x,y) ones(size(x));
% fun_D = inline('ones(size(x))','x','y');
% fun_v = @(x,y) [zeros(size(x)), zeros(size(x))];     % convection term
fun_r = @(x,y) zeros(size(x));                       % reaction term
% fun_r = inline('zeros(size(x))','x','y');
fun_f = @(x,y) 2*pi*pi*sin(pi*x).*sin(pi*y);
% fun_f = inline('2*pi*pi*sin(pi*x).*sin(pi*y)','x','y');
fun_u = @(x,y) sin(pi*x).*sin(pi*y);        %  Analyt
% fun_u = inline('sin(pi*x).*sin(pi*y)','x','y');% solution
fun_g = @(x,y) sin(pi*x).*sin(pi*y);        %  Dirichlet boundary
% fun_g = inline('sin(pi*x).*sin(pi*y)','x','y');

%%% This is another case
% epsilon = 0.05;  
% % u = inline('atan2(y.^2-2*x+1/2,epsilon)','x','y','epsilon');  % or
% fun_u = @(x,y) atan2(y.^2-2*x+1/2,epsilon);
% fun_f = @(x,y) -8*epsilon*(-20*y.*y+24*x-7+4*epsilon*epsilon-12*y.^4+16*y.*y.*x+ ...
%         16*x.*x)./(4*epsilon*epsilon+4*y.^4-16*y.*y.*x+4*y.*y+16*x.*x-8*x+1).^2;
% fun_g = @(x,y) fun_u(x,y);
% fun_a = @(x,y) ones(size(x));
% fun_c = @(x,y) zeros(size(x));
[A,b] = fem_elliptic(dof_map_elem, n_dof, ...
                     Nodes, Elements, ... 
                     fun_D, fun_r, fun_f,  ...                      
                     quad_weight, quad_pnt, ...
                     bas_val, bas_grad_xi, bas_grad_eta);  

%% Apply Dirichlet boundary conditions
bnd_dof = unique([dof_map_edge(DirichletSides,1), dof_map_edge(DirichletSides,2)]);
bnd_val = fun_g(dof_xyz(bnd_dof,1), dof_xyz(bnd_dof,2));
A(bnd_dof,:) = 0;
A(:,bnd_dof) = 0;
A(bnd_dof,bnd_dof) = speye(length(bnd_dof));
%% the same functional with following for - cycle
% for k = 1:length(bnd_dof)
%     A(bnd_dof(k), bnd_dof(k)) = 1;
% end
b(bnd_dof) = bnd_val;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% and solve it
uh = A \ b;  % bicgstab, cgs

% % 7. postprossing: general purpose subroutine
% disp_d = FE_Order; 
% tri_temp = msh_template_tri(disp_d);
% [I,J,K] = indices(disp_d);
% 
% 
% u = fun_u(Nodes(:,1),Nodes(:,2));
% trisurf(Elements,Nodes(:,1),Nodes(:,2),abs(uh-u));

%% by default the first nnode dof is for the nodes
trisurf(Elements(:,1:3),Nodes(:,1),Nodes(:,2),uh(1:size(Nodes,1)));
