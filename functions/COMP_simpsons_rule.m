function [q]=COMP_simpsons_rule(f,a,b,n)
% HELP. COMPOSITE SIMPSON's method. Some numerical calculations
%       and analysis exercises of Numeric Integration. 
%
%       f function is given in terms of a symbolic variable x and as an
%       inline function. E.g., f=inline('x^2+2*x-2'). Also, if the function
%       f is trigonometric function, the 4th argument can be entered as
%       'trigonom' or just 'trig' or 1. X is expected to be in degrees for
%       trigonometric function evaluations. The number of steps NSTEPS has 
%       to be even.
%       upl and lowl are upper and lower limits. NB: A sequential order of 
%       limits is unnecessary to follow, 'if' conditions will take care of
%       lower and upper limits accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Written by Sulaymon L. ESHKABILOV, Ph.D
%                         October, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if numel(f)>1 % If the input provided is a vector
    n=numel(f)-1; h=(b-a)/n;
    I= h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end))+f(end));
else % If the input provided is an anonymous functio
    h=(b-a)/n; xi=a:h:b;xi=xi';
    I= h/3*(f(xi(1))+2*sum(f(xi(3:2:end-2)))+4*sum(f(xi(2:2:end)))+f(xi(end)));
end
q=I;
