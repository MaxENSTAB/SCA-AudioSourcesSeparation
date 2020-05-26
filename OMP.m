function [A,histo] = OMP(D,X,L)
%
%=============================================
% 'Orthogonal matching pursuit' coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
%
% INPUTS :
% D - the dictionary
% X - the signals to represent
% L - the maximal number of coefficient for representation
%     of each signal.
%
% OUTPUTS :
% A - sparse coefficient matrix.
% histo - histogram of the atoms used in the sparse repr.
%=============================================

[n,P]=size(X);
[n,K]=size(D);
A=zeros(K,P);
histo=zeros(1,K);

% Orthogonal matching pursuit
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
    indx=zeros(L,1);
    for j=1:1:L,
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
    end;
    A(indx,k)=a;

end;

% Histogram of atom indices
if (nargout==2)
    for k=1:1:P,
        for i=1:K
            if (A(i,k)~=0)
                histo(i)=histo(i)+1;
            end
        end
    end
end

