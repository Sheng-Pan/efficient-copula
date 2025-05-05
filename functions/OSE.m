function q = OSE(z0,K,paraes0,R1,R2,c,num,comint)
%z0 is the copula family, choices are 'C''Clayton', F'Frank', t't', G'Gumber'
% K is the step number
%paraes0 is the initial estimator
%R1 R2 are the rank of data
%c  and num , comintare tuning parameters in algorithm
%c taylor expansion point, num is the order of taytor expansion, comint is
%the COMPOSITE SIMPSON's integral points number
%e.g copulaeff3(f1,f,PML,R1,R2,0,12,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 f=@(u,v,s)copulapdf(z0,[u*ones(length(v),1),v],s);
     f1=@(u,v,s)copulapdf(z0,[u,v],s);iter=0;
   while(iter<K)
              ab1= copulaeff3(f1,f,paraes0,R1,R2,c,num,comint);
              paraes0 = paraes0 + mean(ab1)/mean(ab1.*ab1);
              iter = iter+1;
   end
q =  paraes0;