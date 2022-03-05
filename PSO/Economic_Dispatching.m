clc;
clear all;
close all;

%%
%problem defenition
global NFE;
NFE=0;
model=CreateModel();

Costfunction=@(p) Cost(p,model);

nvar=model.N;
varsize=[1,nvar];
varmin=model.pmin;
varmax=model.pmax;

%%
%PSO parameters

maxit=1000;

npop=80;


w=1;
wdamp=0.99;
c1=2;
c2=2;

velmax=0.1*(varmax-varmin);
velmin=-velmax;


%%
%Initialization

empty_particle.position=[];
empty_particle.cost=[];
empty_particle.velocity=[];
empty_particle.best.position=[];
empty_particle.best.cost=[];

p=repmat(empty_particle,npop,1);

globalbest.cost=inf;

for i=1:npop
    
   p(i).position=CreateRandomSolution(model);
   
   p(i).velocity=zeros(varsize);
   
   p(i).cost=Costfunction(p(i).position);
   
   p(i).best.position=p(i).position;
   p(i).best.cost=p(i).cost;
   NFE=NFE+1;
   if p(i).best.cost<globalbest.cost
       %globalbest.cost=p(i).best.cost;
       %globalbest.position=p(i).best.cost;
       globalbest=p(i).best;
   end
   
end

bestcost=zeros(maxit,1);
nfe=zeros(maxit,1);

%%
%PSO main loop

for it=1:maxit
    for i=1:npop
        %update velocity
        p(i).velocity=w.*(p(i).velocity)...
            +c1.*rand(varsize).*(p(i).best.position-p(i).position)...
            +c2.*rand(varsize).*(globalbest.position-p(i).position);
        
        %Apply velocity limits
        
        p(i).velocity=min(max(p(i).velocity,velmin),velmax);
        
        %update position
        p(i).position=p(i).position+p(i).velocity;
        
        %Velocity mirror effect
        
        isoutside=(p(i).position<varmin | p(i).position>varmax);
        p(i).position(isoutside)=-p(i).position(isoutside);
        
        %Apply position limits
        p(i).position=min(max(p(i).position,varmin),varmax);
        
        
        
        %evaluation
        p(i).cost=Costfunction(p(i).position);
        NFE=NFE+1;
        
        
        %update personal best
        if p(i).cost<p(i).best.cost
            p(i).best.cost=p(i).cost;
            p(i).best.position=p(i).position;
        
        
            %update global best
            if p(i).best.cost<globalbest.cost
                globalbest=p(i).best;
            end
        end
        
            
    end
    bestcost(it)=globalbest.cost;
    
    nfe(it)=NFE;
    
    disp(['Iteration ' num2str(it) 'NFE = ' num2str(NFE) ': Best Cost = ' num2str(bestcost(it))]);
    
    w=w*wdamp;
end
%%
%Result

figure;
plot(nfe,bestcost,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost');



%%
%Function

function model=CreateModel()
    
    pmin=[808 853 463 857 716 448 539 673];
    pmax=[2166 2172 1526 2177 2166 1788 2041 1513];
    a0=[7067 9488 8882 9702 8213 5175 9161 9577];
    a1=[9 9 9 7 8 6 9 5];
    a2=[-0.1092 -0.1194 -0.2647 -0.2390 -0.1634 -0.2900 -0.1069 -0.1877]*1.0e-4;
    
    N=numel(pmin);
    
    pL=10000;
    
    model.N=N;
    model.pmin=pmin;
    model.pmax=pmax;
    model.a0=a0;
    model.a1=a1;
    model.a2=a2;
    model.pL=pL;
 
end

function p=CreateRandomSolution(model)

    pmin=model.pmin;
    pmax=model.pmax;
    
    p=unifrnd(pmin,pmax);
    
end

function  z=Cost(p,model)
    global NFE;
    
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    z=sum(model.a0+model.a1.*(p)+model.a2.*(p.^2));
    
    %violation
    
    v=abs((sum(p)/model.pL)-1);
    
    beta=10; 
    
    z=z*(1+beta*v);
    
end