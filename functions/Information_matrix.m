function q = Information_matrix(z0,alpha,nu)
%z0 is the copula family, choices are 'C''Clayton', F'Frank', t't', G'Gumber'
%alpha is the copula parameter, nu is the t copula freedom 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-8;
 %% %t efficiency bound 
if z0 == 't'
K=1;N4=10000;f1=@(u,v,s)copulapdf(z0,[u,v],s,nu);
     U1 = copularnd(z0,alpha,nu,N4);
    p00= copulaslcoet(f1,U1,K,alpha,tol);
 q =  copulaslt(f1,p00, alpha,K,U1,tol);
%% %Achimedean efficiency bound 
elseif z0 == 'F'
z0 = 'Frank';
 N4=10000;K=1;f1=@(u,v,s)copulapdf(z0,[u,v],s);
    U1 = copularnd(z0,alpha,N4);
    p00=copulaslcoe(f1,U1,K,alpha,tol);
    q = copulasl(f1,p00, alpha,K,U1,tol);
     while  isnan(q)==1
               U1 = copularnd(z0,alpha,N4);
               p00=copulaslcoe(f1,U1,K,alpha,tol);
               q = copulasl(f1,p00, alpha,K,U1,tol);
      end
elseif z0 == 'C'%clayton
z0 = 'Clayton';
N4=10000;K=2;f1=@(u,v,s)copulapdf(z0,[u,v],s);
    U1 = copularnd(z0,alpha,N4);
    p00=copulaslcoe(f1,U1,K,alpha,tol);
    q = copulasl(f1,p00, alpha,K,U1,tol);
elseif z0 == 'G'%clayton
z0 = 'Gumbel';
N4=10000;K=2;f1=@(u,v,s)copulapdf(z0,[u,v],s);
    U1 = copularnd(z0,alpha,N4);
    p00=copulaslcoe(f1,U1,K,alpha,tol);
    q = copulasl(f1,p00, alpha,K,U1,tol);
 end  

function [q] = copulasl(f,b,s,K,U1,tol)
 u=U1(:,1);v=U1(:,2);%b1=b(1:K);b2=b((K+1):2*K);
  Is=(log(f(u,v,s+tol))-log(f(u,v,s)))/tol;
Iu=(log(f(u+tol,v,s))-log(f(u,v,s)))/tol;
Iv=(log(f(u,v+tol,s))-log(f(u,v,s)))/tol;
aa=Is-spine2(b,u,K)-Iu.* intspine2(b,u,K)-spine2(b,u,K)-Iu.* intspine2(b,u,K);
q=mean(power(aa,2));
function [q1] = copulaslt(f,b,s,K,U1,tol)
 u=U1(:,1);v=U1(:,2);
  Is=(log(f(u,v,s+tol))-log(f(u,v,s)))/tol;
  Iu=(log(f(u+tol,v,s))-log(f(u,v,s)))/tol;
Iv=(log(f(u,v+tol,s))-log(f(u,v,s)))/tol;
aa=Is-spine2(b,v,K)-Iv.* intspine2(b,v,K);
 q1=mean(power(aa,2));
function [a] = spine2(b,u,K)%(x,y,sigma)

a=0;
for k =1:K
    a = a+ b(k).*2^0.5*cos(k*pi*u);
end
function [a] = intspine2(b,u,K)%(x,y,sigma)

a=0;
for k =1:K
    a = a+ b(k).*2^0.5*sin(k*pi*u)/k*pi;
end
function [q] = copulaslcoet(f,U1,K,s,tol)%(x,y,sigma)
 u=U1(:,1);v=U1(:,2);
  Is=(log(f(u,v,s+tol))-log(f(u,v,s)))/tol;
Iv=(log(f(u,v+tol,s))-log(f(u,v,s)))/tol;
w=[];
for k =1:K
    w = [w,2^0.5*cos(k*pi*v)+Iv.*2^0.5.*sin(k*pi*v)/k*pi];
end
q=(w'*w)\w'*Is;
function [q] = copulaslcoe(f,U1,K,s,tol)%(x,y,sigma)
 u=U1(:,1);v=U1(:,2);
  Is=(log(f(u,v,s+tol))-log(f(u,v,s)))/tol;
Iv=(log(f(u,v+tol,s))-log(f(u,v,s)))/tol;
w=[];
for k =1:K
    w = [w,2^0.5*cos(k*pi*v)+Iv.*2^0.5.*sin(k*pi*v)/k*pi];
end
q=(w'*w)\w'*Is;