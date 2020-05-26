function [A,histo] = OMPerr(D,X,resnorm)

%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  resnorm - the norm of the final residue
%
% output arguments: A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);
[n,K]=size(D);
A=zeros(size(D,2),size(X,2));
histo=zeros(K,1);

% Orthogonal Matching Pursuit
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
	indx = [];
	a = [];
	currResNorm2 = sum(residual.^2);
    j=0;
    while currResNorm2>resnorm^2,
		j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
		currResNorm2 = sum(abs(residual).^2);
    end;
   if ~isempty(indx)
       A(indx,k)=a;
   end
end

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

