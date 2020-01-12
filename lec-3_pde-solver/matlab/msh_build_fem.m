function [Nodes, Edges, Elements] = msh_build_fem(V,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [T, E, ET, TE] = msh_build_fem(V,T)
% generate edge information for the simplest mesh V and T, secondly find 
% the neighbouring informations for simple mesh composed of V and T.
% By the way, update the vertices' index into counter clock-wise.
%
% Actually, it is only executed in the preprocessing, so the efficiency
% is not that important, abd we can improve it a little bit later. 
%
% July 14th, 2012.
%
%  Author: Dr. Xian-Liang Hu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the areas of all triangles.
x21 = V(T(:,2),1) - V(T(:,1),1);
x31 = V(T(:,3),1) - V(T(:,1),1);
y21 = V(T(:,2),2) - V(T(:,1),2);
y31 = V(T(:,3),2) - V(T(:,1),2);
areas = 0.5*(x21.*y31 - x31.*y21);
% promising the index of vertices is anti-clock
pos = find(areas<0); % find the anti-direction triangles
tmp = T(pos,2); % and change their position;
T(pos,2) = T(pos,3);
T(pos,3) = tmp;
%areas(pos) = -1*areas(pos);
%% build all edges
E = get_edge_list(T);
%% pre allocate mem space
ne = size(E,1); ET = zeros(ne, 2);
nt = size(T,1); TE = zeros(nt, 3);
%% This is such a slower version, please make it faster by vectorizing the
%% code or c/fortran interface
for elem = 1:nt
    for j = 1:3        
        % now the j'th edge of element elem
        k1 = mod(j ,3) + 1;
        k2 = mod(k1,3) + 1;
        v1 =  T(elem, k1);
        v2 =  T(elem, k2);
        
        % find the eg's global index
        eg = find((v1 == E(:,1))&(v2 == E(:,2)));  %% !!! very bad idea !!!
        if(length(eg) == 1)
            ET(eg, 1) = elem;  % the index in E and T are same
        else
            eg = find((v1 == E(:,2))&(v2 == E(:,1)));
            if(length(eg) == 1)
                ET(eg, 2) = elem;  % the index in E and T are different
            else
                fprintf('There must be something wrong in the mesh!\n');
                return;
            end
        end
        % set the eg to be the j'th edge of triangle elem
        TE(elem,j) = eg;
        
    end
end
% postpross the edge at the boundary
E_filter = ET(:,1) == 0;
ET(E_filter,1) = ET(E_filter,2);
ET(E_filter,2) = 0;
tmp = E(E_filter,1);
E(E_filter,1) = E(E_filter,2);
E(E_filter,2) = tmp;
%% write the output
% The last column of Nodes is the boundary type
nv = size(V,1);
%% [x,y,z, flag];
Nodes = [V zeros(nv,1) zeros(nv,1)];
%% [v1,v2,nei1,nei2, son, daughter, flag];
Edges = [E, ET, -1*ones(ne,2), zeros(ne,1)];
%% [v1,v2,v3, e1, e2, e3, spring, summer, autumn, winter, flag];
Elements = [T, TE, -1*ones(nt,4), zeros(nt,1)];
end
 
%%%%%%%%%%%%%
% E = get_edge_list(T)
function E = get_edge_list(T)
    np = max(max(T));
    ip1 = T(:,1); ip2 = T(:,2); ip3 = T(:,3);
    A = sparse(ip1, ip2, -1, np, np);
    A = A + sparse(ip2, ip3, -1, np, np);
    A = A + sparse(ip3, ip1, -1, np, np); 
    AA = tril(A + A');   % treat the boundary edge, which is not symmetric
    [i,j] = find(AA); 
    E = [i j];
end
	
% function idx = sort_edge(V,T,E,ET,edges,k)
% 
% % find the first ring of vertex k
% nedges = size(edges,1);
% vtx = zeros(ndeges,1);
% for eg = 1:nedges
%     edge = edges(eg);  % global index of the edge
%     vtx(eg) = E(edge,1);
%     if(vtx(eg) == k)
%         vtx(eg) = E(edge,2);
%     end
% end
% 
% % check if it is a boundary vertex
% bnds = find(ET(edges,2)==0);
% if(isempty(bnds))
%     
% else
%     
% end
% 
% 
% % and then sort them in anti-clock turn
% idx = 1:nedges;  
% % idx = 1:2*nedges-1;
% % idx(nedges) = 1;
% % for k = 2:nedges  % sort the other edges according to ET
% %     
% % end
% end