function [Pbest,Cbest,pg] =TLSPSO(f,a,b,D,MaxDT,N,varargin)%◊Ó÷’À„∑®TLS-PSO
%Three-learning strategy particle swarm algorithm,TLSPSO
%_________________________________________________________________________%
%  TLSPSO source codes demo 1.0    %
%  [u0,s,Gbest] =TLSPSO(f,a,b,D,MaxDT,N,varargin)                                    %
% Input:
% f; objective function
% a: the lower boundary vector of the search space
% b: the upper boundary vector of the search sapce
% D: the dimension of the search sapce
% MaxDT: the maximum iteration number
% N: the population size
% Output:
% u0: the best function value obtained by TLSPSO
% Gbest: the global best solution obtained by TLSPSO
% s: the best value record with each iteration
% Developed in MATLAB R2014b                                             %
% Programmer: Xinming Zhang                                   %
%         e-Mail: xinmingzhang@126.com                                    %
%   Main paper:  Xinming Zhang, Qiuying Lin. Three-learning strategy      %
%   particle swarm algorithm. Information Sciences, 2022, 593: 289-313£¨
%   DOI£∫10.1016/j.ins. 2022.01.075                   %

      
c2=1.49445;            %
vmax=max(b-a)*0.1;%
vmin=-vmax;        %
va=vmin*ones(N,D);
vb=vmax*ones(N,D);
%initinize the velocities and positions of the population
v=rand(N,D)*(vmax-vmin)+vmin;
aa=repmat(a,N,1);
bb=repmat(b,N,1);
x=aa+(bb-aa).*rand(N,D);
%-------------
p=feval(f,x',varargin{:});
y=x;
temp=x;
[Pbest,index]=min(p);
pg=x(index,:);%Pg is the global best solution

Cbest=zeros(1,MaxDT);

%------Loop------------
for t=1:MaxDT
    wt=1-(t-1)/MaxDT;     
    w=2-2*t/MaxDT;   
    if t<MaxDT/2
        num=0;
    else
        num=1;
    end
      
    [p,ind]=sort(p);    
    y=y(ind,:);
    v=v(ind,:);
    x=x(ind,:); 
    temp=temp(ind,:);
    for i=num+1:N
        if (i==N)&&(t>=MaxDT*0.9)
            difw=(pg-x(i,:));%Worst-best example learning strategy
            temp(i,:)=c2*rand(1,D).*difw;%
        else            
            [difm,rnum]=MPEL(i,y,x);%Med-point-example learning strategy            
            difr=RL(rnum,i,y,x);%Random learning strategy
            temp(i,:)=(2-w)*rand(1,D).*difm+w*rand(1,D).*difr;
            v(i,:)=wt*v(i,:)+temp(i,:);
        end
    end            
    v=simplebounds(v,va,vb);
    x(1:N,:)=x(1:N,:)+v(1:N,:);
    x=ControlBoundD(x,aa,bb);
    fit=feval(f,x',varargin{:});
    for i=1:N
        value=fit(i);
        if  value< p(i)
            p(i)=value;
            y(i,:)=x(i,:);
            if value < Pbest
                pg=y(i,:);
                Pbest=value;
            end
        end
    end
    Cbest(t)=Pbest;    
end



