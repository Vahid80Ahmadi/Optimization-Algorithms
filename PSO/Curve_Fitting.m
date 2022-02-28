clc;
close all;
clear all;
%%
%Problem  Parameters
global NFE;
NFE=0;


model.x=[0 1 2 4];
model.y=[1 0 2 2];
model.fhat=@(x,a) a(1)+a(2)*x+a(3)*x^2+a(4)*x^3;

cost=@(a) Cost(a,model);

nvar=4;
varsize=[1,nvar];
varmin=-10;
varmax= 10;

%%
%PSO Parameters

iter=600;
npop=100;

w=1;
wp=0.99;
c1=2;
c2=2;
mu=0.05;

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
        
        newsol.position=Mutate(globalbest.position,mu,varmax,varmin);
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
    newsol.position=Mutate(globalbest.position,mu,varmax,varmin);
    [newsol.cost,newsol.sol]=cost(newsol.position);    
    if newsol.cost<=globalbest.cost
       globalbest=newsol;           
    end
       
    w=w*wp;
    
    bestcost(it)=globalbest.cost;
    nfe(it)=NFE;
   
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost = ' num2str(bestcost(it))]);
    
    figure(1);
    Myplot(globalbest.position,model);
end

%%
%Function
function [z,sol]=Cost(a,model)
    
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    x=model.x;
    y=model.y;
    fhat=model.fhat;
    
    yhat=zeros(size(x));    
    for i=1:numel(x)
       yhat(i)=fhat(x(i),a); 
    end
    
    e=y-yhat;
    
    z=sum(e.^2);
    
    sol.yhat=yhat;
    sol.e=e;
end

function Myplot(a,model)
    
    x=model.x;
    y=model.y;
    fhat=model.fhat;
    
    xmin=min(x);
    xmax=max(x);
    dx=xmax-xmin;
    xmin=xmin-0.1*dx;
    xmax=xmax+0.1*dx;
    
    
    xx=linspace(xmin,xmax,100);
    yy=zeros(size(xx));
    for i=1:numel(xx) 
       yy(i)=fhat(xx(i),a); 
    end
    
    plot(x,y,'ro','MarkerFaceColor','y','MarkerSize',10);
    hold on;
    plot(xx,yy,'b','LineWidth',2);
    legend('Main Data','Fitted Curve');
    xlabel('x');
    ylabel('y');
    hold off;
    grid on;
end

function xnew=Mutate(x,mu,varmax,varmin)
    
    nvar=numel(x);
    
    nmu=ceil(mu*nvar);
    
    j=randsample(nvar,nmu);
    
    sigma=0.1*(varmax-varmin);
    
    xnew=x;
    xnew(j)=x(j)+sigma*randn(size(j));
    
    xnew=min(max(xnew,varmin),varmax);
    
end
