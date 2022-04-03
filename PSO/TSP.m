clc;
clear all;
close all;
%%
%Problem  Parameters
global NFE;
NFE=0;


model=CreateModel();

cost=@(x) Cost(x,model);

nvar=model.N;
varsize=[1,nvar];
varmin=0;
varmax=1;

%%
%PSO Parameters

iter=600;
npop=100;

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
        
        newsol.position=Mutate(p(i).position);
        [newsol.cost,newsol.sol]=cost(newsol.position);
        
        if newsol.cost<=p(i).cost
           p(i).position=newsol.position;
           p(i).cost=newsol.cost;
           p(i).sol=newsol.sol;            
        end
          
        if p(i).cost<p(i).best.cost
            
            p(i).best.cost=p(i).cost;
            p(i).best.position=p(i).position;
            p(i).best.sol=p(i).sol;
            
            if p(i).best.cost<globalbest.cost
                
                globalbest=p(i).best;
            end
        end
    end        
    newsol.position=Mutate(globalbest.position);
    [newsol.cost,newsol.sol]=cost(newsol.position);    
    if newsol.cost<=globalbest.cost
       globalbest=newsol;           
    end
       
    w=w*wp;
    
    bestcost(it)=globalbest.cost;
    nfe(it)=NFE;
    
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost = ' num2str(bestcost(it))]);
    
    figure(1);
    GraphPlot(globalbest.sol.Tour,model);
    grid on;
end

%%
%Functions
function model=CreateModel()
    
    x=[10 28 55 96 97 16 98 96 49 81 15 43 92 80 96 66 4 85 94 68];
    y=[76 75 40 66 18 71 4 28 5 10 83 70 32 96 4 44 39 77 80 19];
    
    N=numel(y);
    
    D=zeros(N,N);
    
    for i=1:N-1
        for j=i:N            
            D(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            D(j,i)=D(i,j);
        end
    end
    
    model.N=N;
    model.x=x;
    model.y=y;
    model.D=D;
    
end

function [z,sol]=Cost(x,model)

    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    L=0;
    N=model.N;
    D=model.D;
    [~,Tour]=sort(x);
    
    for k=1:N
        
       i=Tour(k); 
        
        if k<N
            j=Tour(k+1);
        else
            j=Tour(1);
        end
        
        L=L+D(i,j);
    end
    z=L;
    
    sol.Tour=Tour;
    sol.L=L;
    
end
function GraphPlot(tour,model)
    
    x=model.x;
    y=model.y;
    tour=[tour tour(1)];
    
    plot(x(tour),y(tour),'bs-','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','y');
end

function xnew=Mutate(x)
    
    [~,tour]=sort(x);
    
    M=randi([1,3]);
    
    switch M
        case 1
            t=mswap(tour);
        case 2
            t=mrev(tour);
        case 3
            t=minser(tour);
    end
    
    xnew=zeros(size(x));
    xnew(t)=x(tour);
end

function t=mswap(tour)
    
    n=numel(tour);
    
    i=randsample(n,2);
    i1=i(1);
    i2=i(2);
    
    t=tour;
    t([i1 i2])=tour([i2 i1]);
    
end

function t=mrev(tour)
    
    n=numel(tour);
    
    i=randsample(n,2);
    i1=min(i);
    i2=max(i);
    
    t=tour;
    t(i1:i2)=tour(i2:-1:i1);

end

function t=minser(tour)
    
    n=numel(tour);
    
    i=randsample(n,2);
    i1=i(1);
    i2=i(2);
    
    if i1>i2       
        t=[tour(1:i2-1) tour(i2+1:i1) tour(i2) tour(i1+1:end)];
    else
        t=[tour(1:i1) tour(i2) tour(i1+1:i2-1) tour(i2+1:end)];
    end

end