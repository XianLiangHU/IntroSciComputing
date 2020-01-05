function [A, b] = fem_elliptic(dof_map, n_dof, ...
                                Nodes, Elements, ...    % mesh
                                fun_D, fun_r, fun_f, ...                % param_fun
                                quad_weight, quad_pnt,   ...            % quadrature
                                bas_val, bas_grad_xi, bas_grad_eta)     % finite element basis                  
         
%% do the elementary calculation with gien mesh and fe and coef functions

nt = size(Elements,1);
n_bas = size(bas_val,2);
% reserve the mem space for equation system
ii = repmat(dof_map, n_bas, 1);  ii = reshape(ii, nt, n_bas*n_bas);
jj = repmat(dof_map, 1, n_bas);  jj = reshape(jj, nt, n_bas*n_bas);
A_nnz = zeros(nt, n_bas*n_bas);
b_nnz = zeros(nt, n_bas);

%% when nonlinear triangle is used, it should be evaled point-wisely
x21 = Nodes(Elements(:,2),1) - Nodes(Elements(:,1),1);
y21 = Nodes(Elements(:,2),2) - Nodes(Elements(:,1),2);
x31 = Nodes(Elements(:,3),1) - Nodes(Elements(:,1),1);
y31 = Nodes(Elements(:,3),2) - Nodes(Elements(:,1),2);
Jacob = x21.*y31 - x31.*y21; 

%% calculate the geometrical factor
dxi_dx  =  y31./Jacob;  deta_dx = -y21./Jacob;
dxi_dy  = -x31./Jacob;  deta_dy =  x21./Jacob;

for k = 1:length(quad_weight)
    qxw = 0.5*quad_weight(k)*Jacob;  % quad weight on triangles
    
    %% eval the gradients at the current quadrature points
    dphi_dx = dxi_dx*bas_grad_xi(k,:) + deta_dx*bas_grad_eta(k,:);
    dphi_dy = dxi_dy*bas_grad_xi(k,:) + deta_dy*bas_grad_eta(k,:);
    
    %% eval the rhs at the current quadrature points
    px = Nodes(Elements(:,1),1) + x21*quad_pnt(k,1) + x31*quad_pnt(k,2);
    py = Nodes(Elements(:,1),2) + y21*quad_pnt(k,1) + y31*quad_pnt(k,2);
    f_val = feval(fun_f, px, py);
    D_val = feval(fun_D, px, py);
    r_val = feval(fun_r, px, py);   
    
    for i = 1:n_bas
        for j = 1:n_bas
            ij = n_bas*(i - 1) + j;
            %% diffusion term:
            A_val = D_val.*(dphi_dx(:,i).*dphi_dx(:,j) + dphi_dy(:,i).*dphi_dy(:,j));
            %% reaction term:
            A_val = A_val + r_val*(bas_val(k,i)*bas_val(k,j));
            %% do the numerical integration:
            A_nnz(:, ij) = A_nnz(:,ij) + qxw.*A_val;   
        end
        %% do intergration for right hand side:
        b_nnz(:,i) = b_nnz(:,i) + qxw.*(f_val*bas_val(k,i));
    end
    
end

% generate the sparse linear system
A = sparse(ii, jj, A_nnz, n_dof, n_dof);
b = full(sparse(dof_map, ones(nt, n_bas), b_nnz, n_dof, 1));

return