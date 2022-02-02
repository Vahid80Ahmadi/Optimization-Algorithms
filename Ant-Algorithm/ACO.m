clc;
clear all; 
close all;
%% Problem Parameters

model=CreateModel();

x=model.x;
y=model.y;
D=model.D;
n=model.n;

cost=@(x) Cost(x,D);

%% ACO Parameters

maxiter=200;

npop=50;

alpha=1;
beta=1;

rho=0.5;

tau0=0.8;

Q=0.6;
%% Setup ACO

tau=tau0*ones(n,n);
eta=1./D;

ant.path=[];
ant.cost=[];
pop=repmat(ant,npop,1);

bestsol.path=[];
bestsol.cost=inf;

bestcost=zeros(maxiter,1);
%% Main Loop

for it=1:maxiter
  
    for k=1:npop
        
        pop(k).path=randi([1 n]);
        
        for l=2:n
            
            i=pop(k).path(end);
            
            p=(tau(i,:).^alpha).*(eta(i,:).^beta);
            
            p(pop(k).path)=0;
            
            p=p/sum(p);
            
            j=Roulett(p);
            
            pop(k).path=[pop(k).path j];
            
        end
        
        pop(k).cost=cost(pop(k).path);
        
        if pop(k).cost<bestsol.cost
            bestsol=pop(k);
        end
        
    end
        
    %Update Phromones

    for k=1:n
        tour=pop(k).path;

        tour=[tour tour(1)];

        for l=1:n
            i=tour(l);
            j=tour(l+1);
            
            tau(i,j)=tau(i,j)+Q/pop(k).cost;
            
        end

    end    
    
    % Evaporation
    
    tau=tau*(1-rho);
    
    bestcost(it)=bestsol.cost;
    
    disp(['iteration: ' num2str(it) ' Cost: ' num2str(bestcost(it))]);
    
    
   % pplloott(bestsol.path,model);
    
end


%% Functions
function model=CreateModel()
    
    x=[66 4 85 94 68 76 75 40 66 18 71 4 28 5 10 83 70 32 96 4];
    y=[44 39 77 80 19 49 45 65 71 76 28 68 66 17 12 50 96 35 59 23];
    
    n=numel(x);
    
    D=zeros(n);
    
    for i=1:n
        for j=i:n
            D(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            D(j,i)=D(i,j);
        end
    end
    
    
    model.D=D;
    model.x=x;
    model.y=y;
    model.n=n;
    
end

function z=Cost(x,D)
    
    z=0;
    x=[x x(1)];
    n=numel(x);
    
    for i=1:n-1
        z=z+D(x(i),x(i+1));
    end

end

function i=Roulett(p)
    
    r=rand;
    c=cumsum(p);
    i=find(r<=c,1,'first');
end


function pplloott(path,model)
    
    figure(1)
    path=[path path(1)];
    
    x=model.x;
    y=model.y;
    
    plot(x,y,'o','MarkerFaceColor','y','MarkerSize',9);
    hold on;
    
    X=x(path);
    Y=y(path);
    
    figure(1)
    plot(X,Y);
        
    hold off;
end
