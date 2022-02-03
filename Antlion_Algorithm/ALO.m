clc;
clear all;
close all;
%% Problem Parameters

cost=@(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmin=-10;
varmax= 10;

lb=varmin*ones(varsize);
ub=varmax*ones(varsize);
%% ALO Parameters

nAnt=50;

nAntLion=round(0.6*40);

maxit=500;

I=2;    % First Value for I

%% Setup ALO

empty_ant.position=[];
empty_ant.cost=[];

empty_antlion.position=[];
empty_antlion.cost=[];

ant=repmat(empty_ant,nAnt,1);
antlion=repmat(empty_antlion,nAntLion,1);

for i=1:nAnt
    
    ant(i).position=unifrnd(varmin,varmax,varsize);
    ant(i).cost=cost(ant(i).position);

end

for i=1:nAntLion
    
    antlion(i).position=unifrnd(varmin,varmax,varsize);
    antlion(i).cost=cost(ant(i).position);

end

% ant=Sorting(ant);
antlion=Sorting(antlion);
bestcosts=zeros(maxit,1);
%% Main Loop

for it=1:maxit
    
   % Defininig w  
   if it>0.1*maxit
       w=2;
   elseif it>0.5*maxit
       w=3;
   elseif it>0.75*maxit
       w=4;
   elseif it>0.9*maxit
       w=5;
   else
       w=6;
   end
   
   % Defining I 
   I=(10^w)*(it/maxit);
   
   % Update Bounds
   c=varmin*ones(varsize)/I;
   d=varmax*ones(varsize)/I;
   
   % Random Walking
   % AntLions
   
   % AntLion Selection Perobability
   AntLionCosts=[antlion.cost];
   AntLionCosts=1./AntLionCosts;
   pAntLion=AntLionCosts/(sum(AntLionCosts));
   
   for i=1:nAnt
        j=Roulett(pAntLion);
        
        if rand<0.5
            lb=lb+antlion(j).position;
        else
            lb=-lb+antlion(j).position;
        end
        
        if rand>=0.5
            ub=ub+antlion(j).position;
        else
            ub=-ub+antlion(j).position;
        end
        
        ant(i).position=((ant(i).position-max(ant(i).position)).*(d-c))/...
            (max(ant(i).position)-min(ant(i).position))+c;
        
        ant(i).position=min(max(ant(i).position,c),d);
        
        ant(i).cost=cost(ant(i).position);
        
        if ant(i).cost<antlion(j).cost
            antlion(j)=ant(i);
        end
        
   end
   
   antlion=Sorting(antlion);
   
   bestsol=antlion(1);
   bestcosts(it)=bestsol.cost;
   
   disp(['Iteration: ' num2str(it) ' Cost: ' num2str(bestcosts(it))]);
   
end

%% Funcitons
function z=Sphere(x)
    z=sum(x.^2);
end

function pop=Sorting(pop)
    [~,so]=sort([pop.cost]);
    pop=pop(so);
end

function i=Roulett(p)
    
    r=rand;
    c=cumsum(p);
    
    i=find(r<=c,1,'first');
end