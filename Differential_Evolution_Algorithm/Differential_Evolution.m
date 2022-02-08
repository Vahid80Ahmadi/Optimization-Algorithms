clc;
clear all;
close all;
%% 
%Problem Parameters

cost=@(x) Sphere(x);

nvar=20;
varsize=[1 nvar];

varmin=-5;
varmax=5;

%%
%DE Parameters

maxit=1000;
npop=50;

beta_min=0.2;
beta_max=0.8;

pcr=0.2;

%%
%Setup DE

empty.position=[];
empty.cost=[];
%empty.diff=[];
bestsol.cost=inf;
bstsol.position=[];
pop=repmat(empty,npop,1);

for i=1:npop
    
    pop(i).position=unifrnd(varmin,varmax,varsize);
    pop(i).cost=cost(pop(i).position);
    
    if pop(i).cost<bestsol.cost
        bestsol=pop(i);
    end
end

bestcost=zeros(maxit,1);

%%
%Main Loop

for it=1:maxit
    
    for j=1:npop
        
        x=pop(j).position;
        A=randperm(npop);
        
        A(A==j)=[];
        a=A(1);
        b=A(2);
        c=A(3);
        
        %Mutation
        MyBeta=unifrnd(beta_min,beta_max,varsize);
        y=pop(a).position+(MyBeta.*(pop(b).position-pop(c).position));
        y=min(max(y,varmin),varmax);
        %Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        
        for i=1:numel(x)           
            if j==j0 || rand<=pcr
                z(i)=y(i);
            else
                z(i)=x(i);
            end 
        end
        
        newsol.position=z;
        newsol.cost=cost(newsol.position);
        
        if newsol.cost<pop(j).cost
           pop(j)=newsol;
           if pop(j).cost<bestsol.cost
               bestsol=pop(j);
           end
        end
        
    end
    
    bestcost(it)=bestsol.cost;
    
    
    disp(['Iteration : ' num2str(it) ' BestCost : ' num2str(bestcost(it))]);
    
end

%%
%Results

figure;
plot(bestcost);

%%
%Functions
function z=Sphere(x)
    z=sum(x.^2);
end