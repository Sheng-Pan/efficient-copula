function q = OSEt(K,paraes0,nu,R1,R2,num)
%z0 is the copula family, choices are 'C''Clayton', F'Frank', t't', G'Gumber'
% K is the step number
%paraes0 is the initial estimator
%nu freedom parameter
%R1 R2 are the rank of data
%c  and num , comintare tuning parameters in algorithm
%c taylor expansion point, num is the order of taytor expansion, comint is
%the COMPOSITE SIMPSON's integral points number
%e.g copulaeff3(f1,f,PML,R1,R2,0,12,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter=0;
   while(iter<K)
              ab1= copulaeff(num,paraes0,nu,R1,R2);
              paraes0 = paraes0 + mean(ab1)/mean(ab1.*ab1);
              paraes0(paraes0 >= 1) = 0.99;
              paraes0(paraes0  <= -1) = -0.99;
              iter = iter+1;
   end
q =  paraes0;
function [q] = copulaeff(n,s,nu,R1,R2)
%calculation method of efficient score of t copula
% n COMP_simpsons_rule number 
% s copula parameter nu freedom parameter
% R1 ,R2 input rank of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=@(u,v,s)copulapdf('t',[u*ones(length(v),1),v],s,nu);
f1=@(u,v,s)copulapdf('t',[u,v],s,nu);
N=999;
tol=1e-6;%u=0.2;v=0.4;s=2;
a=1e-3;b=1-1e-3;h=(b-a)/n; xi=a:h:b;alpha1=zeros(n+1,1);gamma1=zeros(n+1,1);ker=zeros(length(xi)+2,length(xi)+2);
r1=zeros(length(R1),1);r2=zeros(length(R1),1);
for k=1:n+1
a1=zeros(length(xi),1)';
a1(3:2:end-2)=4*kernelcopula(f1,xi(k).*ones(length(xi(3:2:end-2)),1),xi(3:2:end-2)',s,tol)';
a1(2:2:end)= 2*kernelcopula(f1,xi(k).*ones(length(xi(2:2:end)),1),xi(2:2:end)',s,tol)';
    ker(k+1,:)=[1,a1,1]; 
   alpha1(k)=     alphacopula(f,xi(k),s,tol,N);
   gamma1(k)=  gamcopula(f,xi(k),s,N,tol);
end
for k=1:length(R1)
      [~,s1] = min(abs(R1(k) - xi));
      [~,s2] = min(abs(R2(k) - xi));
       r1(k)=s1;
       r2(k)=s2;
end
alpha1=diag([0;alpha1;0]);n0=n+3;gamma1=[0;gamma1;0];
b1 = n^2*ones(n0-1,1); b2 =  - 2*n^2*ones(n0,1);
C = diag(b1,-1) + diag(b2,0) + diag(b1,+1);C(1,:) = 0;C(n0,:) = 0;C(1,1) = 1;C(n0,n0) = 1;
H=(C-alpha1+ker)\(-gamma1);
h=[0;diff(H)/n];
Is= (log(f1(R1,R2,s+tol))-log(f1(R1,R2,s)))/tol;
Iu=(log(f1(R1+tol,R2,s))-log(f1(R1,R2,s)))/tol;Iv=(log(f1(R1,R2+tol,s))-log(f1(R1,R2,s)))/tol;
q=Is-h(r1)-h(r2)-Iu.*H(r1)-Iv.*H(r2);

function [q] = alphacopula(f,u,s,tol,N)

q=zeros(length(u),1);
 for k=1:length(u)
  fun=@(v)alphacopula1(f,u(k),v,s,tol);
q(k)=COMP_simpsons_rule(fun,0.0001, 0.9999,N);
end

function [q] = alphacopula1(f,u,v,s,tol)
 Iu= (log(f(u+tol,v,s))-log(f(u,v,s)))/tol;
q=power(Iu,2).*f(u,v,s);

function [q] = gamcopula(f,u,s,N,tol)
q=zeros(length(u),1);
 for k=1:length(u)
   fun=@(v)gamcopula1(f,u(k),v,s,tol);
q(k)=COMP_simpsons_rule(fun,0.0001, 0.9999,N);
end



function [q] = gamcopula1(f,u,v,s,tol)%(x,y,sigma)
 Is= (log(f(u,v,s+tol))-log(f(u,v,s)))/tol;
 Iu= (log(f(u+tol,v,s))-log(f(u,v,s)))/tol;
q=Is.*Iu.*f(u,v,s);

function [q] = kernelcopula(f,u,v,s,tol)%(x,y,sigma)
fun=@(v) (log(f(u+tol,v,s))-log(f(u,v,s)))/tol;
q=((fun(v+tol)-fun(v))/tol).*f(u,v,s);





