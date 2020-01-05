function [phi, dphidxi, dphideta] = fe_lagrange(x, y, FE_Order)
% this funtion return the value of shape functions of degree d at points v
%  d ------ degree , must be positive integer
%  x and y is the corresponding value of L1 and L2
%  the master element         < (0,0),(1,0),(0,1)> 
% the derivtive must be comput with relate to the jacobi of element

J  = x;
K = y;
I = 1 - J - K;

%% This is a lazy version
phi = vdm23(FE_Order, I, J, K);

desc = desc_pattern(FE_Order + 2);

Id = diag(ones((FE_Order + 1)*(FE_Order + 2)/2,1));
xi = FE_Order*de_cast_step(Id, FE_Order,-1,1,0, desc);  % direction at v1->v2
eta = FE_Order*de_cast_step(Id,FE_Order,-1,0,1,desc);  % direction v1->v3

Mat = vdm23(FE_Order-1,I,J,K);  % low order vdm

dphidxi = Mat*xi;  %get the deriv B-form value
dphideta = Mat*eta;

return



%%
function Bout = de_cast_step(Bin,d,lam1,lam2,lam3,desc_pattern)
% This function does one step of the de Casteljau algorithm
% Note: This function assumes that the input arguments are all column vectors
% global desc_pattern;
m_in = size(Bin,1);
% d = degree(m_in);

%% below is the alternatable algorithms
m_out = m_in-d-1;
n = length(lam1);
indx1 = desc_pattern(1:m_out,1);  % always 1:m_out, but we place it here for simplicity
indx2 = desc_pattern(1:m_out,2);
indx3 = desc_pattern(1:m_out,3);
if size(Bin,2) == 1
   Bout = Bin(indx1)*lam1' + Bin(indx2)*lam2' + Bin(indx3)*lam3';
else
   Bout = Bin(indx1,:)*spdiags(lam1,0,n,n) + Bin(indx2,:)*spdiags(lam2,0,n,n) + ...
      Bin(indx3,:)*spdiags(lam3,0,n,n);
end;

% % the former algorithm is a little faster
% m_rows = d*(d+1)/2;
% m_cols = (d+1)*(d+2)/2;
% I = (1:m_rows)';
% Id = ones(m_rows,1);
% desc = sparse(I,desc_pattern(I,1),lam1*Id,m_rows,m_cols) + ...
%         sparse(I,desc_pattern(I,2),lam2*Id,m_rows,m_cols) + ...
%         sparse(I,desc_pattern(I,3),lam3*Id,m_rows,m_cols);
% Bout = desc*Bin;
return

%%
function pattern = desc_pattern(d)
% this function may be execute only once in whole comput session.
% much like the function asce_pattern.
m = (d+1)*(d+2)/2;
pattern = zeros(m,3);
begin = 1;
for j = 0:d
    idx = (begin:(begin+j))';
    pattern(idx,1) = idx;
    pattern(idx,2) = idx + j + 1;
    pattern(idx,3) = idx + j + 2;
    begin = begin + j + 1;
end
return

%%
function Mat = vdm23(d,b1,b2,b3)
% comput the Bform of degreed d at points b1,b2,b3;
m = (d+1)*(d+2)/2;
plot_m = length(b1);
[I,J,K] = indices(d);
IM = diag(I)*ones(m,plot_m);
JM = diag(J)*ones(m,plot_m);
KM = diag(K)*ones(m,plot_m);

plot_IM = diag(b1)*ones(plot_m,m);
plot_JM = diag(b2)*ones(plot_m,m);
plot_KM = diag(b3)*ones(plot_m,m);
Mat = (plot_IM).^(IM').*(plot_JM).^(JM').*(plot_KM).^(KM');
IF = gamma(I+1);
JF = gamma(J+1);
KF = gamma(K+1);
A = factorial(d)*ones(plot_m,m)*diag(1./(IF.*JF.*KF));
Mat = A.*Mat;
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