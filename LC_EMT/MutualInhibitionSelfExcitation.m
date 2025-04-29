function out = MutualInhibitionSelfExcitation(t, x, P)
%This function is a tristable system representing a two genes which inhibit
%the other's expression, while promoting their own. 

%Note that everything in this function is vectorized. 
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


[num,m]=size(x);
H=zeros(num,m,num); 
for i=1:16
    for j=1:16 
        if IM(i,j)==1
            H(i,:,j)=Ha(x(i,:),S(i,j),AB(i,j),Hilln(i,j));
        elseif IM(i,j)==-1
            H(i,:,j)=Hr(x(i,:),S(i,j),AB(i,j),Hilln(i,j));
        end
    end
end
k=ones(1,16);
% k(1)=P(1);k(4)=P(2);k(6)=P(3);k(8)=P(4);k(12)=P(5);k(13)=P(6);k(14)=P(7);k(15)=P(8);
k(1:9)=P(1:9);k(11:15)=P(10:14);
for i=1:16
    F(i,:)=sum(H(:,:,i))-k(i)*x(i,:);
end

%Baseline values for P: [.5 .5 .5 .5]; 

out = F;
end
function H=Hr(X,S,b,n)
H=b.*S.^n./(X.^n+S.^n);
end
function H=Ha(X,S,a,n)
H=a.*X.^n./(S.^n+X.^n);
end
   