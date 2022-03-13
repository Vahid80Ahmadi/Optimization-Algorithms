clc;
clear all;
close all;
%%
%Problem  Parameters
global NFE;
NFE=0;

model=CreateModel();

cost=@(x) BCost(x,model);

nvar=model.N;
varsize=[1,nvar];
varmin=zeros(varsize);
varmax=ones(varsize);

%%
%PSO Parameters

iter=600;
npop=200;

w=1;
wp=0.99;
c1=2;
c2=2;

velmax=0.1*(varmax-varmin);
velmin=-velmax;

%%
%Setup PSO
particle.position=[];
particle.velocity=[];
particle.cost=[];
particle.sol=[];
particle.best.position=[];
particle.best.cost=[];
particle.best.sol=[];

globalbest.position=[];
globalbest.cost=inf;

p=repmat(particle,npop,1);

for i=1:npop
   
    p(i).position=unifrnd(varmin,varmax,varsize);
    p(i).velocity=zeros(varsize);
    [p(i).cost,p(i).sol]=cost(p(i).position);
    p(i).best.sol=p(i).sol;    
    p(i).best.position=p(i).position;
    p(i).best.cost=p(i).cost;
   
    if p(i).best.cost<globalbest.cost
        globalbest=p(i).best;
    end
end


%%
%PSO Main Loop
bestcost=zeros(iter,1);
nfe=zeros(iter,1);

for it=1:iter
    
    for i=1:npop
        p(i).velocity=w*p(i).velocity...
            +c1*rand(varsize).*(p(i).best.position-p(i).position) ...
            +c2*rand(varsize).*(globalbest.position-p(i).position);
        
        p(i).velocity=min(max(p(i).velocity,velmin),velmax);
        
        p(i).position=p(i).position+p(i).velocity;

        isout=(p(i).position<varmin | p(i).position>varmax);
        p(i).velocity(isout)=-p(i).velocity(isout);
        
        p(i).position=min(max(p(i).position,varmin),varmax);
        
        [p(i).cost,p(i).sol]=cost(p(i).position);
        
%         newsol.position=Mutate(globalbest.position,mu,varmax,varmin);
%         [newsol.cost,newsol.sol]=cost(newsol.position);
%         
%         if newsol.cost<=p(i).cost
%            p(i).position=newsol.position;
%            p(i).cost=newsol.cost;
%            p(i).sol=newsol.sol;            
%         end
%           
        if p(i).cost<p(i).best.cost
            
            p(i).best.cost=p(i).cost;
            p(i).best.position=p(i).position;
            p(i).best.sol=p(i).sol;
            
            if p(i).best.cost<globalbest.cost
                
                globalbest=p(i).best;
            end
        end
    end        
%     newsol.position=Mutate(globalbest.position,mu,varmax,varmin);
%     [newsol.cost,newsol.sol]=cost(newsol.position);    
%     if newsol.cost<=globalbest.cost
%        globalbest=newsol;           
%     end
       
    w=w*wp;
    
    bestcost(it)=globalbest.cost;
    nfe(it)=NFE;
   
    
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost = ' num2str(bestcost(it))]);
    

end
    figure;
    plot(nfe,bestcost,'LineWidth',2);
%%
%Function
function model=CreateModel()
    
    model.mi=[5 1 2 3 7 1 1 2 6 6 6 4];
    model.val=[9 2 5 3 4 5 9 1 19 19 10 10];
    model.wi=[23 46 25 14 41 25 19 26 13 15 48 49];
    model.w=1000;
    model.N=numel(model.wi);
    
end

function [z,sol]=FCost(x,model)
        
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    varmax=model.mi;
    val=model.val;
    wi=model.wi;
    w=model.w;
    
    z1=sum(val.*(varmax-x));
 
    v=max((sum(wi.*x)/w)-1,0);
    beta=1000;

    z=z1+beta*v;
    
    sol.v=v;
    sol.z1=z1;
    sol.z=z;
    sol.weight=sum(wi.*x);
end

function [z,sol]=BCost(x,model)
        
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;

    nvar=model.N;
    varsize=[1,nvar];
    varmin=zeros(varsize);
    varmax=model.mi;
    
    k=min(floor(varmin+(varmax-varmin+1).*x),varmax);
    
    val=model.val;
    wi=model.wi;
    w=model.w;
    
    z1=sum(val.*(varmax-k));
 
    v=max((sum(wi.*k)/w)-1,0);
    beta=1000;

    z=z1+beta*v;
    
    sol.position=k;
    sol.val=val;
    sol.v=v;
    sol.z1=z1;
    sol.z=z;
    sol.weight=sum(wi.*k);
    
end