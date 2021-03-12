% Chuong trinh Lanczos tinh cho ma tran A doi xung
% Kiem tra OK
clc, clear all
b = [1;2;2]
A = [-1,2,-3;-5,3,2;1,5,5];
A = A'*A
q(:,1) = b/norm(b)
beta(1) = 0;
m = 3;
for j = 1:m
    if j == 1
        j
        v(:,j)=A*q(:,j)
        anpha(j) = v(:,j)'*q(:,j)
        v(:,j) = v(:,j) - anpha(j)*q(:,j)
        beta(j+1) = norm(v(:,j))
        q(:,j+1) = v(:,j)/beta(j+1)
    elseif j > 1
        j
        v(:,j) = A*q(:,j) - beta(j)*q(:,j-1)
        anpha(j) = v(:,j)'*q(:,j)
        v(:,j) = v(:,j) - anpha(j)*q(:,j)
        beta(j+1) = norm(v(:,j))
        if beta(j+1)<10^(-20)
            break
        end
        q(:,j+1) = v(:,j)/beta(j+1)
    end
end
T = diag(anpha);
T = T+ diag(beta(2:length(beta)-1),1);
T = T+ diag(beta(2:length(beta)-1),-1);
T
be = norm(b);
e1 = zeros(length(T),1);
e1(1,1) = 1;
y = inv(T)*be*e1
[r,c] = size(q);
q = q(1:r,1:c-1)
x = q*y
