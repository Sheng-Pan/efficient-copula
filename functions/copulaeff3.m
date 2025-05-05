function [q] = copulaeff3(f1,f,s,R1,R2,c,Xnumber,intenumber)
%efficient score calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = Xnumber;a=1e-3;b=1-1e-3;aa = 0:1/N:1;%aa(1)=a;aa(end)=b;
N1=intenumber;tol=1e-8;h=(b-a)/N1; xi=a:h:b;
H1 = [];alpha1=zeros(N+1,1);gamma1=zeros(N+1,1);intkern=zeros(N+1,N+1);
for m =1:(N+1)
    H1 = [H1,power(aa-c,m-1)'];%x-c (x-c)^2...
    alpha1(m)=     alphacopula(f,aa(m),s,tol,10);
    gamma1(m)=  gamcopula(f,aa(m),s,N,tol);
    for j = 1:(N+1)
        fun = @(x)kernelcopula(f1,x.*ones(length(xi(3:2:end-2)),1),xi(3:2:end-2)',s,tol).*...
                power(x-c,j-1)/ factorial(j-1);
        a1 = 4*fun(aa(m))';
        a2 = 2*fun(aa(m))';
        intkern(m,j) = a1(1) + a2(1);%simpson integral
    end
end
xx = [ones(N+1,1),repmat(power(cumprod(1:N),-1),N+1,1)];%1 1/2! 1/3!...
H1 = H1.*xx;%x-c/1 (x-c)^2/2!...
Alpha1 = repmat(alpha1,1,N+1).*H1;
Alpha1(1,:)=0;intkern(1,:)=0;gamma1(1)=0;%H(0)=0
Alpha1((N+1),:)=0;intkern((N+1),:)=0;gamma1((N+1))=0;%H(1)=0
H2=[zeros(N+1,1),zeros(N+1,1),H1(:,1:(N-1))];H2(1,:)=H1(1,:);H2(N+1,:)=H1(N+1,:);
H0=(H2-Alpha1+intkern)\(-gamma1);
H = @(x)[ones(length(x),1),cumprod(repmat(x-c,1,N),2).*repmat(power(cumprod(1:N),-1),length(x),1)]*...
    H0;
h = @(x)[ones(length(x),1),cumprod(repmat(x-c,1,N-1),2).*repmat(power(cumprod(1:N-1),-1),length(x),1)]*...
    H0(2:end);
Is= (log(f1(R1,R2,s+tol))-log(f1(R1,R2,s)))/tol;
Iu=(log(f1(R1+tol,R2,s))-log(f1(R1,R2,s)))/tol;Iv=(log(f1(R1,R2+tol,s))-log(f1(R1,R2,s)))/tol;
q=Is-h(R1)-h(R2)-Iu.*H(R1)-Iv.*H(R2);
function [q] = kernelcopula(f,u,v,s,tol)%(x,y,sigma)
fun=@(v) (log(f(u+tol,v,s))-log(f(u,v,s)))/tol;
q=((fun(v+tol)-fun(v))/tol).*f(u,v,s);
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



