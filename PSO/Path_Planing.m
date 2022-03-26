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

iter=500;
npop=100;

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
        
        %Y 
        p(i).velocity.y=w*p(i).velocity.y...
            +c1*rand(varsize).*(p(i).best.position.y-p(i).position.y) ...
            +c2*rand(varsize).*(globalbest.position.y-p(i).position.y);
        
        p(i).velocity.y=min(max(p(i).velocity.y,velmin.y),velmax.y);
        
        p(i).position.y=p(i).position.y+p(i).velocity.y;

        isout=(p(i).position.y<varmin.y | p(i).position.y>varmax.y);
        p(i).velocity.y(isout)=-p(i).velocity.y(isout);
        
        p(i).position.y=min(max(p(i).position.y,varmin.y),varmax.y);
        
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
    Picture(globalbest.sol,model);
end

%%
%Functions
function model=CreateModel()
    
    xs=0;
    ys=0;
    
    xt=4;
    yt=6;
    
    xc=3;
    yc=4;
    r=2;
    
    xmax= 10;
    xmin=-10;
    ymax= 10;
    ymin=-10;
    
    n=5;
    
    model.xs=xs;
    model.ys=ys;
    model.xt=xt;
    model.yt=yt;
    model.xc=xc;
    model.yc=yc;
    model.r=r;
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
    
    xs=model.xs;
    ys=model.ys;
    xt=model.xt;
    yt=model.yt;
    
    XS=[xs x xt];
    YS=[ys y yt];
    
    k=numel(XS);
    
    TS=linspace(0,1,k);
    
    tt=linspace(0,1,100);
    xx=spline(TS,XS,tt);
    yy=spline(TS,YS,tt);
    
    xc=model.xc;
    yc=model.yc;
   
    
    dx=diff(xx);
    dy=diff(yy);
    
    L=sum(sqrt(dx.^2+dy.^2));
    d=sqrt((xx-xc).^2+(yy-yc).^2);
    
    sol2.xx=xx;
    sol2.yy=yy;
    sol2.XS=XS;
    sol2.YS=YS;
    sol2.L=L;
    sol2.d=d;
end

function Picture(sol,model)

    xs=model.xs;
    ys=model.ys;
    xt=model.xt;
    yt=model.yt;
    xc=model.xc;
    yc=model.yc;
    r=model.r;
    
    XS=sol.XS;
    YS=sol.YS;
    xx=sol.xx;
    yy=sol.yy;
    
    theta=linspace(0,2*pi,100);
    fill(xc+r*cos(theta),yc+r*sin(theta),'r');
    hold on;
    plot(xx,yy,'b','LineWidth',2);
    plot(XS,YS,'ro');
    plot(xs,ys,'bs','MarkerSize',12,'MarkerFaceColor','y');
    plot(xt,yt,'bs','MarkerSize',12,'MarkerFaceColor','y');
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
    
    z=sol2.L;
    
    d=sol2.d;
    r=model.r;
    
    v=max(1-(d/r),0);
    vt=mean(v);
    beta=10;
    
    z=z*(1+beta*vt);
    
    sol.z=z;
    sol.v=vt;
    sol.XS=sol2.XS;
    sol.YS=sol2.YS;
    sol.xx=sol2.xx;
    sol.yy=sol2.yy;
    
end