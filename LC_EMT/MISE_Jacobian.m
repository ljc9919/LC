function out = MISE_Jacobian(t, xin, P)
if ndims(xin)<3
    in = zeros(1, size(xin, 1), size(xin, 2)); 
    in(1, :, :) = xin; 
    xin=in; 
end
tsteps = size(xin, 3); %If there is a constant or a zero in the jacobian...
format long
 
% x1 = xin(:, 1, :); 
% x2 = xin(:, 2, :); 
%  
% n = 4; k1 = 1; k2 = 1;  a1 = 1; a2 = 1; b1 = .1; b2 = .1;
%  
% out = [x1.^(-1).*(P(1).^n+x1.^n).^(-2).*(a1.*n.*P(1).^n.*x1.^n+(-1).*k1.*x1.*(P(1).^n+x1.^n).^2),...
%     (-1).*b1.*n.*P(2).^n.*x2.^((-1)+n).*(P(2).^n+x2.^n).^(-2);...
%     (-1).*b2.*n.*P(4).^n.*x1.^((-1)+n).*(P(4).^n+x1.^n).^(-2),...
%     x2.^(-1).*(P(3).^n+x2.^n).^(-2).*(a2.*n.*P(3).^n.*x2.^n+(-1).*k2.*x2.*(P(3).^n+x2.^n).^2)];
 
IM=load ('matrix.txt');
N = size(IM,1); %The dimension of the system.
Hilln=zeros(N,N); %Hill functiom
S=zeros(N,N); % threshold
AB=zeros(N,N); %the scale factors for the activation and inhibition
for i=1:N
    for j=1:N
        Hilln(i,j)=4;
        if IM(i,j)==1 && i~=j
            AB(i,j) = 2; %a
            if i>10 && j<11
            S(i,j) = .22; %sb
            else
            S(i,j) = 6; %sa
            end
        elseif IM(i,j)==1 && i==j
            AB(i,j) = 10; %aa
            if i>10 && j<11
            S(i,j) = .22; %sb
            else
            S(i,j) = 6; %sa
            end
        elseif IM(i,j)==-1
            AB(i,j) = 2; %b
            if i>10 && j<11
            S(i,j) = .22; %sb
            else
            S(i,j) = 6; %sa
            end
        end
    end
end


num=size(xin,2);

Jacobi_multi = zeros(num,num,tsteps);
k=ones(1,16);
% k(1)=P(1);k(4)=P(2);k(6)=P(3);k(8)=P(4);k(12)=P(5);k(13)=P(6);k(14)=P(7);k(15)=P(8);
k(1:9)=P(1:9);k(11:15)=P(10:14);
for i = 1:tsteps
    jacobi_tensor = Jacobi_add(xin(:,:,i),k,S,IM,Hilln,AB);
    Jacobi_multi(:,:,i) = jacobi_tensor;
end

out = Jacobi_multi;
end

function jacobi_matrix = Jacobi_add(x,k,S,IM,Hilln,AB)
num=size(IM,1);
jacobi_matrix = zeros(num,num);
for i = 1:num
    for j = 1:num
        if IM(i,j)==0
            J_ij=0;
        end
        if IM(i,j)==1
            J_ij =  AB(i,j)*deri_Ha(x(i),S(i,j),Hilln(i,j));
        end
        if IM(i,j)==-1
            J_ij =  AB(i,j)*deri_Hi(x(i),S(i,j),Hilln(i,j));
        end
        jacobi_matrix(i,j) = J_ij;
    end
end

jacobi_matrix = jacobi_matrix - k.*eye(num);

end

function deri_value = deri_Ha(x,s,n)
deri_value = s.^n.*n.*x.^(n-1)./(s.^n+x.^n).^2;
end

function deri_value = deri_Hi(x,s,n)
deri_value = -s.^n.*n.*x.^(n-1)./(s.^n+x.^n).^2;
end

