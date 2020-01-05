function [dof_map_elem, dof_map_side, dof_xyz, dof_patten_side, dof_xieta_elem, n_dof] = dof_lagrange(Nodes, Edges, Elements, FE_Order)
% function [dof_elem, dof_edge, dof_xyz, dof_patten_side, n_dof] = dof_lagrange(V,T,E, FE_Order)
%    generate necessary informations about degrees of freedom for mesh
%    (V,T,E) and FE_Order order lagrange element
%
%% sort dofs in c^0 style between elements
% but in a hierarchy style within each element
% that is to say, each column of dof(:,i) means dofs for element i in order
%       1
%       2 3
%       4 5 6
%       7 8 9 10
%       ... ... ...  

    [dof_map_elem, dof_map_side, n_dof] = build_dof_map(Nodes,Edges, Elements,FE_Order);
%     dof_edge = build_dof_edge(dof_elem, Elements, Edges);
    dof_xyz  = build_dof_xyz (Nodes, Elements, dof_map_elem, n_dof, FE_Order);
    
    dof_patten_side = element_dof_pattern(FE_Order) ;    
    [I, J, K] = indices(FE_Order);
    dof_xieta_elem = [J, K]/FE_Order;
return

%%  distribute dof for high order Lagrange/Bernstein element
function [dof_map_elem, dof_map_side, n_dof] = build_dof_map(Nodes, Edges, Elements, FE_Order)
%% reserve space for output
n_dof = 0;

%% 1. build dofs for vertices
n_vtx = size(Nodes,1);
dof_map_vtx = n_dof + (1:n_vtx)';
n_dof = n_vtx;

%% 2. sord dofs to edges
n_edge = size(Edges,1);
dof_map_side = zeros(n_edge, FE_Order+1); %% same with Edges in linear FE case
dof_map_side(:, 1) = dof_map_vtx(Edges(:, 1));
line_dof = (FE_Order - 1)*((1:n_edge)' - 1);
for k = 2:FE_Order
   dof_map_side(:, k) = (n_dof + k - 1) + line_dof;
end
dof_map_side(:, FE_Order+1) = dof_map_vtx(Edges(:, 2));
n_dof = n_dof + (FE_Order - 1)*n_edge;

% 3. build dof_map for each elements
n_elem = size(Elements,1);
m = (FE_Order + 1)*(FE_Order + 2)/2;  %% sane with Elements in linear FE case
dof_map_elem = zeros(n_elem, m);

% 3.1 from nodes
dof_map_elem(:, 1) = dof_map_vtx(Elements(:, 1));
dof_map_elem(:, m - FE_Order) = dof_map_vtx(Elements(:, 2));
dof_map_elem(:, m) = dof_map_vtx(Elements(:, 3));

% 3.2 from edges, only for degree > 1 FE
if (FE_Order > 1)
 
    pat = element_dof_pattern(FE_Order);
    all_edge = (1:n_edge)';
    
    %
    for k = 1:3
        %% % should not be 0, which means outside of domain
        T_left = Edges(:,3);  
        E_left = all_edge(:);
        eg_k = (Elements(T_left, 3+k) == E_left);  % pick out the coincide one
        dof_map_elem(T_left(eg_k), pat((2:FE_Order)',k)) = ...
            dof_map_side(E_left(eg_k), (2:FE_Order)');
        
        %% there may be boundary element
        % the dof index in right neighbour must be the inverse of the index
        % in this edge
        eg_filter = Edges(:, 4) ~= 0;
        T_right = Edges(eg_filter, 4)';
        E_right = all_edge(eg_filter);
        eg_k = (Elements(T_right, 3+k) == E_right);  % pick out the coincide one
        dof_map_elem(T_right(eg_k), pat((2:FE_Order)', k)) = ...
            dof_map_side(E_right(eg_k), (FE_Order:-1:2)');
    end
    
end

% 3.3 inside element
if (FE_Order > 2) 
    % find out the inner index for high order FE,
    % logically relative to element <V1,V2,V3>
    % it is enssential to C^r (r >= 1) smooth condition
    m_inner = (FE_Order - 2)*(FE_Order - 1)/2;  
    patt = zeros(m_inner,1); 
    patt(1) = 5;
    pos = 1;
    for i = 1:(FE_Order - 3)
        for j = 0:i
            patt(pos + j + 1) = patt(pos) + j + 3;
        end
        pos = pos + i + 1;
    end
    % 
    dof_map_elem(:,patt) = reshape((n_dof + 1) : (n_dof + n_elem*m_inner), n_elem, m_inner);
    n_dof = n_dof + n_elem*m_inner;
end

return

%%%%  a local index, be careful to chang it!
function pat = element_dof_pattern(FE_Order) 
    pat = zeros(FE_Order + 1,3); %
    m = (FE_Order + 1)*(FE_Order + 2)/2;
    pat(:,1) = ((m - FE_Order) : m)';  % last line
    pat(1,3) = 1;     % first one of line 3
    pat(FE_Order + 1,2) = 1;  % last one of line 2
    for i = 1:FE_Order
       pat(i+1,3) = pat(i,3) + i;  % sort from low to high
       pat(FE_Order + 1 - i,2) = pat(FE_Order + 2 - i,2) + (i + 1);  % from high to low
    end
return

%% %%%%%%%%%%%%%%%%%%%%%%%%
function dof_xyz = build_dof_xyz(Nodes, Elements, dof_map_elem, n_dof, FE_Order)   
DOW = 2;
dof_xyz = zeros(n_dof,DOW);

[I, J, K] = indices(FE_Order);
I = I/FE_Order; J = J/FE_Order; K = K/FE_Order;

m = (FE_Order + 1)*(FE_Order + 2)/2;
n_elem = size(Elements, 1);
for k = 1:DOW
    idx = reshape(dof_map_elem, m*n_elem, 1);
    dof_xyz(idx, k) = reshape(Nodes(Elements(:, 1), k)*I' + ...
                              Nodes(Elements(:, 2), k)*J' + ...
                              Nodes(Elements(:, 3), k)*K', m*n_elem, 1);
end
return

function [I,J,K] = indices(d)
%        [I,J,K] = indices(d)
m = (d+1)*(d+2)/2;
I = zeros(m,1);
J = I;
K = I;
idx = 1;
for i = d:(-1):0
    for j = (d-i):(-1):0
        I(idx) = i;
        J(idx) = j;
        K(idx) = d - i - j;
        idx = idx + 1;
    end
end
return