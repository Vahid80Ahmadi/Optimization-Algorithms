clc;
close all;
clear all;
%%
%Peroblems Parameters

cost=@(x) MOP2(x);

nvar=3;
varsize=[nvar 1];

varmax= 4;
varmin=-4;

nobj=numel(cost(unifrnd(varmin,varmax,varsize)));

%%
%MOEAD Parameters

maxit=200;
npop=50;
nArchive=50;

T=max(ceil(0.2*npop),2);

params.gamma=0.5;
params.varmin=varmin;
params.varmax=varmax;

%%
%Setup MOEAD

%Initialize Sub Problem
sp=CreateSubProblem(nobj,npop,T);

%Initialize Goal Point
z=zeros(nobj,1);

%Initialize Population
empty_individual.position=[];
empty_individual.cost=[];
empty_individual.g=[];
empty_individual.isdominated=[];

pop=repmat(empty_individual,npop,1);

for i=1:npop
   pop(i).position=unifrnd(varmin,varmax,varsize); 
   pop(i).cost=cost(pop(i).position);
   z=min(z,pop(i).cost);
end

for i=1:npop
   pop(i).g=DecomposedCost(pop(i),z,sp(i).lambda);
end

pop=DetermineDomination(pop);
EP=pop(~[pop.isdominated]);
%%
%Main Loop

for it=1:maxit
    
   for i=1:npop 
    
       %Reproduction (Crossover)
       k=randsample(T,2);
       i1=sp(i).neighbors(k(1));
       i2=sp(i).neighbors(k(2));
       
       p1=pop(i1);
       p2=pop(i2);
       
       y=empty_individual;
       y.position=Crossover(p1.position,p2.position,params);
       y.cost=cost(y.position);
       
       z=min(z,y.cost);
       
       for j=sp(i).neighbors  
          y.g=DecomposedCost(y,z,sp(j).lambda); 
          if y.g<=pop(j).g
             pop(j)=y;
          end
       end

   end

   pop=DetermineDomination(pop);

   ndpop=pop(~[pop.isdominated]);

   EP=[EP
       ndpop];
   
   if numel(EP)>nArchive
        Extra=numel(EP)-nArchive;
        ToBeDeleted=randsample(numel(EP),Extra);
        EP(ToBeDeleted)=[];
   end
   
   EP=DetermineDomination(EP);
   EP=EP(~[EP.isdominated]);

   EPC=[EP.cost];

   figure(1);
   plot(EPC(1,:),EPC(2,:),'ro');
   
   
end
%%
%Functions

function z=MOP2(x)
    
    n=numel(x);
    
    z=[1-exp(-sum((x-1/sqrt(n)).^2))
       1-exp(-sum((x+1/sqrt(n)).^2))];
end

function sp=CreateSubProblem(nobj,npop,T)

    empty_sp.lambda=[];
    empty_sp.neighbors=[];

    sp=repmat(empty_sp,npop,1);

    for i=1:npop
        lambda=rand(nobj,1);
        lambda=lambda/norm(lambda);
        sp(i).lambda=lambda; 
    end

    LAMBDA=[sp.lambda]';

    D=pdist2(LAMBDA,LAMBDA);

    for i=1:npop
        [~,so]=sort(D(i,:));
        sp(i).neighbors=so(1:T);  
    end

end

function g=DecomposedCost(x,z,lambda)    
    g=max(lambda.*abs(x.cost-z));
end

function y=Crossover(x1,x2,params)
    
    gamma=params.gamma;
    varmin=params.varmin;
    varmax=params.varmax;
    
    alpha=unifrnd(-gamma,1+gamma,size(x1));
    
    
    y=alpha.*x1+(1-alpha).*x2;
    
    y=min(max(y,varmin),varmax);
end

function b=Dominate(x,y)
    b=(all(x.cost<=y.cost) && any(x.cost<y.cost));
end

function pop=DetermineDomination(pop)
   
    npop=numel(pop);
    
    for i=1:npop
       pop(i).isdominated=false;
    end
    
    for i=1:npop
        for j=i+1:npop
            if Dominate(pop(i),pop(j))
                pop(j).isdominated=true;
            elseif Dominate(pop(j),pop(i))
                pop(i).isdominated=true;
            end
        end
    end
    
    
end