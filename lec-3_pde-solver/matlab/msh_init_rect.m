function [V, T] = msh_init_rect(xx,yy,mesh_type)
% function [V, T] = msh_init_rect(xx,yy,mesh_type)
% generate triangulations for rectangular domain
% mesh_type could be any one of following to indicate the directions
%           'type3A'  'type3B', 'type4A', 'type4B'
% xx,yy:  one dimensional nodes lists along x- and y- axis,
%            the nonuniformity of xx and yy yields different mesh in 2D
%
% Example:
%  1. Chebyshev nodes:
%        n=10; k = 0:n;
%        x = 0.5+0.5*cos(k*pi/n);
%        y = 0.5+0.5*cos(k*pi/n);
%        [V, T] = init_mesh_rectangle(x,y,'type3C');
%
%
%
%


nelx = length(xx) - 1;
nely = length(yy) - 1;

%%%%%%%%%%the node coordinate in Cartesian mesh%%%%%%
%
% 11 - 12- 13- 14- 15
%  |   |   |   |   | 
%  6 - 7 - 8 - 9 - 10
%  |   |   |   |   |
%  1 - 2 - 3 - 4 - 5 
%
V = zeros((nelx+1)*(nely+1),2);%the number of grid points
for  i=1:nely+1
   V((i-1)*(nelx+1)+1:i*(nelx+1),1) = xx;
   V((i-1)*(nelx+1)+1:i*(nelx+1),2) = yy(i);
end

%%%%%%% Fill triangle elememt %%%%%%%

if(strcmp(mesh_type,'type3A'))
%  ---------------------
%  |13\14|15\16| 17\18 | 
%  ---------------------
%  | 7\8 | 9\10| 11\12 | 
%  ---------------------
%  | 1\2 | 3\4 |  5\6  |
%  ---------------------
    T = zeros(2*nelx*nely, 3); %Index of nodes constructing each element
    idx = 1;
    for j = 1:nely
        for i = 1:nelx
            T(idx,:)    = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1];
            T(idx+1,:) =  [j*(nelx+1)+i+1 j*(nelx+1)+i (j-1)*(nelx+1)+i];
            idx = idx + 2;
        end
    end

elseif(strcmp(mesh_type,'type3B'))
    T = zeros(2*nelx*nely, 3); %Index of nodes constructing each element
    idx = 1;
    for j = 1:nely
        for i = 1:nelx
            T(idx,:)    = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i];
            T(idx+1,:) =  [j*(nelx+1)+i+1 j*(nelx+1)+i (j-1)*(nelx+1)+i+1];
            idx = idx + 2;
        end
    end
    
elseif(strcmp(mesh_type,'type3C'))
    %    ---------------  
    %     |    / |     /|
    %    |   /  |   /  |
    %    |  /   | /    |
    %    --------------- load here
    %    | \    | \    |
    %    |  \   |   \  |
    %    |    \ |     \|
    %    --------------- 
    T = zeros(2*nelx*nely, 3); %Index of nodes constructing each element
    idx = 1;
    for j = 1:nely/2
        for i = 1:nelx       
            T(idx,:)    = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i];
            T(idx+1,:) =  [j*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1];      
            idx = idx + 2;
        end
    end

    for j = nely/2+1:nely
        for i = 1:nelx
            T(idx,:)    = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1];
            T(idx+1,:) =  [j*(nelx+1)+i+1 j*(nelx+1)+i (j-1)*(nelx+1)+i];
            idx = idx + 2;
        end
    end


elseif(strcmp(mesh_type,'type3D'))
    %    ---------------  
    %     |    / |     /|
    %    |   /  |   /  |
    %    |  /   | /    |
    %    --------------- load here
    %    | \    | \    |
    %    |  \   |   \  |
    %    |    \ |     \|
    %    --------------- 
    T = zeros(2*nelx*nely, 3); %Index of nodes constructing each element
    idx = 1;
    for j = 1:2:nely
        for i = 1:nelx       
            T(idx,:)    = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i];
            T(idx+1,:) =  [j*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1];      
            idx = idx + 2;
        end
    end

    for j = 2:2:nely
        for i = 1:nelx
            T(idx,:)    = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1];
            T(idx+1,:) =  [j*(nelx+1)+i+1 j*(nelx+1)+i (j-1)*(nelx+1)+i];
            idx = idx + 2;
        end
    end
    
elseif(strcmp(mesh_type,'type4A'))
    %    ---------------
    %    | \    |     /|
    %    |  \ 7 | 6 /  |
    %    | 8  \ | /  5 |
    %    ---------------
    %    |    / | \ 4  |
    %    |1  /  |   \  |
    %    |  / 2 | 3   \|
    %    ---------------   
    T = zeros(2*nelx*nely, 3); %Index of nodes constructing each element
    idx = 1;
    for j = 1:2:nely
        for i = 1:2:nelx       
            T(idx,:) = [(j-1)*(nelx+1)+i j*(nelx+1)+i+1 j*(nelx+1)+i];
            T(idx+1,:) =  [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1];
            T(idx+2,:) = [(j-1)*(nelx+1)+i+1 (j-1)*(nelx+1)+i+2 j*(nelx+1)+i+1];
            T(idx+3,:) = [(j-1)*(nelx+1)+i+2 j*(nelx+1)+i+2 j*(nelx+1)+i+1];
            T(idx+4,:) = [j*(nelx+1)+i+1 j*(nelx+1)+i+2 (j+1)*(nelx+1)+i+2];
            T(idx+5,:) = [j*(nelx+1)+i+1 (j+1)*(nelx+1)+i+2 (j+1)*(nelx+1)+i+1];
            T(idx+6,:) = [j*(nelx+1)+i+1 (j+1)*(nelx+1)+i+1 (j+1)*(nelx+1)+i];
            T(idx+7,:) = [j*(nelx+1)+i+1 (j+1)*(nelx+1)+i j*(nelx+1)+i];
            idx = idx+8;
        end
    end

elseif(strcmp(mesh_type,'type4B'))
    % for centers of the grid points
    p = zeros(nelx*nely,2);
    for i = 1:nely   
        p((i-1)*nelx+1:i*nelx,1) = 0.5*(xx(1:end-1)+xx(2:end));
        p((i-1)*nelx+1:i*nelx,2) = 0.5*(yy(i) + yy(i+1));
    end
    V = [V; p];
    
    T = zeros(4*nelx*nely, 3); %Index of nodes constructing each element
    idx = 1;
    %    ---------------
    %    | \   3     / |
    %    |  \      /   |
    %    |    \  /     |
    %    |  4   +   2  |
    %    |    /   \    |
    %    |   /  1   \  |
    %    | /          \|
    %    ---------------  
    for j = 1:nely
        for i = 1:nelx       
            T(idx,:) = [(j-1)*(nelx+1)+i (j-1)*(nelx+1)+i+1 (nely+1)*(nelx+1)+(j-1)*nelx+i];
            T(idx+1,:) = [(j-1)*(nelx+1)+i+1 j*(nelx+1)+i+1 (nely+1)*(nelx+1)+(j-1)*nelx+i];
            T(idx+2,:) = [j*(nelx+1)+i+1 j*(nelx+1)+i (nely+1)*(nelx+1)+(j-1)*nelx+i];
            T(idx+3,:) = [j*(nelx+1)+i (j-1)*(nelx+1)+i (nely+1)*(nelx+1)+(j-1)*nelx+i];
            idx = idx + 4;
        end
    end


else
    
    
end