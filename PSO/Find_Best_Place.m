clc;
clear all;
close all;
%%
%Problem  Parameters
global NFE;
NFE=0;

model=CreateModel();

cost=@(x) Cost(x,model);

nvar=model.n;
varsize=[1,nvar];
varmin.x=model.xmin;
varmax.x=model.xmax;

varmin.y=model.ymin;
varmax.y=model.ymax;

%%
%PSO Parameters

iter=100;
npop=50;

w=1;
wp=0.99;
c1=2;
c2=2;

velmax.x=0.1*(varmax.x-varmin.x);
velmin.x=-velmax.x;

velmax.y=0.1*(varmax.y-varmin.y);
velmin.y=-velmax.x;

%%
%Setup PSO
list.x=zeros(1,npop);
list.y=zeros(1,npop);

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
   
    p(i).position=CreateRandomSolution(model);
    list.x(i)=p(i).position.x;
    list.y(i)=p(i).position.y;
    p(i).velocity.x=zeros(varsize);
    p(i).velocity.y=zeros(varsize);
    [p(i).cost,p(i).sol]=cost(p(i).position);
    p(i).best.sol=p(i).sol;    
    p(i).best.position=p(i).position;
    p(i).best.cost=p(i).cost;
    
    if p(i).best.cost<globalbest.cost
        globalbest=p(i).best;
    end
end
    figure(1);
    Picture(list,model);

%%
%PSO Main Loop
bestcost=zeros(iter,1);
nfe=zeros(iter,1);

for it=1:iter
    
    for i=1:npop
        
        %X
        p(i).velocity.x=w*p(i).velocity.x...
            +c1*rand(varsize).*(p(i).best.position.x-p(i).position.x) ...
            +c2*rand(varsize).*(globalbest.position.x-p(i).position.x);
        
        p(i).velocity.x=min(max(p(i).velocity.x,velmin.x),velmax.x);
        
        p(i).position.x=p(i).position.x+p(i).velocity.x;

        isout=(p(i).position.x<varmin.x | p(i).position.x>varmax.x);
        p(i).velocity.x(isout)=-p(i).velocity.x(isout);
        
        p(i).position.x=min(max(p(i).position.x,varmin.x),varmax.x);
        list.x(i)=p(i).position.x;
        %Y 
        p(i).velocity.y=w*p(i).velocity.y...
            +c1*rand(varsize).*(p(i).best.position.y-p(i).position.y) ...
            +c2*rand(varsize).*(globalbest.position.y-p(i).position.y);
        
        p(i).velocity.y=min(max(p(i).velocity.y,velmin.y),velmax.y);
        
        p(i).position.y=p(i).position.y+p(i).velocity.y;

        isout=(p(i).position.y<varmin.y | p(i).position.y>varmax.y);
        p(i).velocity.y(isout)=-p(i).velocity.y(isout);
        
        p(i).position.y=min(max(p(i).position.y,varmin.y),varmax.y);
        list.y(i)=p(i).position.y;
        
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
    
    figure(1);
    Picture(list,model);
end

%%
%Functions
function model=CreateModel()
    
    xp=[-5 4 7 -2 2 4];
    yp=[3 -3 7 0 4 0];
    
    xmax= 10;
    xmin=-10;
    ymax= 10;
    ymin=-10;
    
    n=1;
    
    model.xp=xp;
    model.yp=yp;
    model.xmin=xmin;
    model.xmax=xmax;
    model.ymax=ymax;
    model.ymin=ymin;
    model.n=n;
end

function sol1=CreateRandomSolution(model)
    
    n=model.n;
    
    xmin=model.xmin;
    xmax=model.xmax;
    ymin=model.ymin;
    ymax=model.ymax;
    
    x=unifrnd(xmin,xmax,[1 n]);
    y=unifrnd(ymin,ymax,[1 n]);
    
    sol1.x=x;
    sol1.y=y;
    
end

function sol2=SecondOrderSolution(sol1,model)
    
    x=sol1.x;
    y=sol1.y;
    xp=model.xp;
    yp=model.yp;
    
    d=sum(sqrt((x-xp).^2+(y-yp).^2));
    
    sol2.d=d;
    sol2.x=x;
    sol2.y=y;
end

function Picture(sol,model)
    
    xp=model.xp;
    yp=model.yp;
    
    x=sol.x;
    y=sol.y;
    
    plot(xp,yp,'bs','MarkerFaceColor','y','MarkerSize',12);
    hold on;
    plot(x,y,'ko');
    hold off;
    grid on;
    axis equal;
    
end 

function [z,sol]=Cost(sol1,model)
    
    global NFE;
    if isempty(NFE)
        NFE=0;
    end 
    NFE=NFE+1;
    
    sol2=SecondOrderSolution(sol1,model);
    
    z=sol2.d;
    
    sol.x=sol2.x;
    sol.y=sol2.y;
    
end