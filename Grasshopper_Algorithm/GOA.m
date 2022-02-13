clc;
clear;
close all;
%% Problem Parameters

cost=@(x) Sphere(x);

nvar=3;

varsize=[1 nvar];

varmin=-10*ones(varsize);
varmax= 10*ones(varsize);

%% GOA Paremeters

maxit=550;

npop=30;

cmin= 1e-5;
cmax= 1;

%% Setup GOA

empty.position=[];
empty.cost=[];

pop=repmat(empty,npop,1);

bestsol=empty;

bestsol.cost=inf;

for i=1:npop
    
    pop(i).position=unifrnd(varmin,varmax);
    pop(i).cost=cost(pop(i).position);
    
    if pop(i).cost < bestsol.cost
        bestsol=pop(i);
    end
    
end

bestcosts=zeros(maxit,1);

%% Main Loop

for it=1:maxit

    c = cmax - it*((cmax-cmin)/maxit);
    
    for i=1:npop
        
        SI = SocialIntraction(pop,npop,varmin,varmax,varsize,nvar,c,i);
        
        pop(i).position = c * SI + bestsol.position;
        
        pop(i).position=min(max(pop(i).position,varmin),varmax);
        
        pop(i).cost=cost(pop(i).position);
        
        if pop(i).cost < bestsol.cost
            bestsol=pop(i);
        end
        
    end
    
    bestcosts(it)=bestsol.cost;
    
    disp(['Iteration : ' num2str(it) ' Cost : ' num2str(bestcosts(it))]);
end

%% Functions

function z=Sphere(x)
    z=sum(x.^2);
end

function SI=SocialIntraction(pop,npop,varmin,varmax,varsize,nvar,c,i)
    
    pos=[pop.position];
    pos=reshape(pos,npop,nvar);
    
    x=pos(i,:);
    pos(i,:)=[];
    
    d=pdist2(x,pos);
    d=2+2*d/max(d);
    
    SI=zeros(varsize);
    
    s=SFunction(d);
    s=s./d;
    
    for j=1:npop-1
        
        SI = SI +0.5*c*((varmin-varmax)*(s(j))).*(x-pos(j,:));
        
    end

    
end

function s=SFunction(d)
    
    f=0.8;
    l=1.3;
    
    s = f*exp(-d/l) - exp(-d);
end