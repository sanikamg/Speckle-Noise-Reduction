function denoised = funETV(noisy,maxiter,lambda,mu,nu)
%% This function solves mixed noise problem using split-Bregman
% Problem formulation is :
% min_(x,s)  |y-x-s|_2+ lambda |s|_1 + mu |Dh*X*D|_1 + mu |Dv*XD|_1
%here y is noisy image and Dh,Dv are horizontal and vertical
% finite difference operators, s is sparse noise, D is 1D total
%variation.

%% Some initializations
[m,n,dim]=size(noisy); mn=m*n;
y=reshape(noisy,mn*dim,1);
b1=zeros(mn*dim,1); x=y;b2=b1; p=b1;q=b1; s=b1;

% Make total variation matrix
Dh=TVmatrix(m,n,'H');
Dv=TVmatrix(m,n,'V');
D=opTV1(dim);D=D';
D1=kron(D',Dh); D2=kron(D',Dv);
%D1=KronProd({Dh,D'}); D2=KronProd({Dv,D'});
%% Main iteration
for i=1:maxiter
    
    %solve subproblem for x
    bigY=(y-s)+nu*D1'*(p-b1)+nu*D2'*(q-b2);
    [x,~]=lsqr(@afun,bigY,1e-15,10,[],[],x);
    
    p=SoftTh(D1*x+b1,mu/nu);
    q=SoftTh(D2*x+b2,mu/nu);
    s=SoftTh(y-x,lambda);
    
    %Update B1,B2 and B3
    
    b1=b1+D1*x-p;
    b2=b2+D2*x-q;
    
    if rem(i,10)==0
        fprintf('%d iterations done of %d \n ',i,maxiter);
    end
end
denoised=reshape(x,m,n,dim);
%% This is a function handle used in LSQR
    function y = afun(x,str)
        tempval= nu*((D1'*(D1*x))+(D2'*(D2*x)))+ x;
        switch str
            case 'transp'
                y = tempval;
            case 'notransp'
                y = tempval;
        end
    end

end
%% Soft Thresholding
function X=SoftTh(B,lambda)

X=sign(B).*max(0,abs(B)-(lambda/2));

end
%%  Total Variation
%This function will generate total variation operator for an image of size
%mxn. Both horizontal and vertical total variation operators can be made
%using this code.
function opD=TVmatrix(m,n,str)

if str=='H' % This will give matrix for Horizontal Gradient
    D = spdiags([-ones(n,1) ones(n,1)],[0 1],n,n);
    D(n,:) = 0;
    D = kron(D,speye(m));
    
elseif str=='V' %This will give matrix for Verticle Gradient
    D = spdiags([-ones(m,1) ones(m,1)],[0 1],m,m);
    D(m,:) = 0;
    D = kron(speye(n),D);
end
opD=D;

end
%% This function will generate total variation operator for 1D signals
function opD=opTV1(m)

%Make two vectors of elements -1 and 1 of lengths (m-1)
B=[ -1*ones(m-1,1),ones(m-1,1)];
%Make sparse diagonal matrix D of size m-1 by m
%with -1's on zeroth diagonal and 1's on 1st diagonal;
D=spdiags(B,[0,1],m-1,m);
%add a last row of all zeros in D
D(m,:)=zeros(1,m);
%Make it as operator
opD=D; %It will convert to operator.
end

