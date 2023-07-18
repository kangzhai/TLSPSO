function TLSPSOmain()
%Three-learning strategy particle swarm algorithm,TLSPSO
%_________________________________________________________________________%
%  TLSPSO source codes demo 1.0    %
%  [u0,s,Gbest] =TLSPSO(f,a,b,D,NEF,N)                                    %
% Input:
% f; objective function
% a: the lower boundary vector of the search space
% b: the upper boundary vector of the search sapce
% D: the dimension of the search sapce
% NEF: the maximum number of evaluation functin
% N: the population size
% Output:
% u0: the best function value obtained by TLSPSO
% Gbest: the global best solution obtained by TLSPSO
% s: the best value with 
%  Developed in MATLAB R2014b                                             %
%  Programmer: Xinming Zhang                                   %
%         e-Mail: xinmingzhang@126.com                                    %
%   Main paper:  Xinming Zhang, Qiuying Lin.      %
%Differential mutation and novel social learning particle optimization    %
%algorithm. Information Sciences, 2022, 593: 289-313£¬                     %
%DOI£º10.1016/j.ins. 2022.01.075                                          %
%                                                                         %
%_________________________________________________________________________%
clc;
clear;
Num=51;
U(1:14)=-1400:100:-100;U(15:28)=100:100:1400;
D=30;
MaxDT=3000;
N=100;
func_num=10;
Vt=zeros(1,Num);
f=str2func('cec13_func');
a=-100*ones(1,D);b=100*ones(1,D);
u1=U(func_num);
CBest=zeros(1,MaxDT);
time=0;
for i=1:Num
    tic;
    [u0,s,Gbest] =TLSPSO(f,a,b,D,MaxDT,N,func_num);
    time=((i-1)*time+toc)/i;
    CBest=CBest+s;
    Vt(i)=u0;
end
Vt=Vt-u1;
MeanValue=mean(Vt);
StdValue=std(Vt);
GoodValue=min(Vt);
BadValue=max(Vt);
plot(Vt)
title([' DSPSO£ºMean=',num2str(MeanValue),'£¬Std=',num2str(StdValue)]);
xlabel(['Vma=',num2str(BadValue),'£¬Vmi=',num2str(GoodValue),'£¬Time=',num2str(time)]);
figure(2)
plot(CBest/Num-u1);
xlabel('Iteration')
ylabel('AverageError')


