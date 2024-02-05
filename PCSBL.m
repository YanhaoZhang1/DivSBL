function x_new=PCSBL(y,A,sigma,beta)

%%
% Input: 
% y: measurements;
% A: sensing matrix;
% sigma: square root of the noise covariance;
% beta: parameter controling the relevance between the elements;
% Output:
% x_new: estimated sparse signal

%%

[m,n]=size(A);
iter=0;
iter_mx=300;
D=eye(n);
sigma2=sigma^2;
alpha_new=ones(n,1);
var_new=inv(A'*A/sigma2+D);
mu_old=ones(n,1);
mu_new=1/sigma2*var_new*A'*y;
a=0.5;
b=1e-10;

while iter<iter_mx& norm(mu_new-mu_old)>1e-6
    iter=iter+1;
    mu_old=mu_new;
    mul=[mu_new(2:n);0];
    mur=[0;mu_new(1:n-1)];
    var=diag(var_new);
    varl=[var(2:n);0];
    varr=[0;var(1:n-1)];
    E=mu_new.^2+beta*mul.^2+beta*mur.^2+var+beta*varl+beta*varr;
    alpha_new=a./(0.5*E+b);
    idx1=find(alpha_new>1e10);
    alpha_new(idx1)=1e10;
    alf=[alpha_new(2:n); 0];                                %   left-shifted version
    arf=[0; alpha_new(1:n-1)];                              %   right-shifted version
    D=diag(alpha_new+beta*alf+beta*arf);
    var_new=inv(A'*A/sigma2+D);
    mu_new=1/sigma2*var_new*A'*y;
end

x_new=mu_new;