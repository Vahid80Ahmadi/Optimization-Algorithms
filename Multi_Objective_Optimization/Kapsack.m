function [z,sol,vol]=Kapsack(x,model)
 
    %n=model.n;
    W=model.W;
    wi=model.wi;
    vi=model.vi;
    Mi=model.Mi;
    
    z1=sum(vi.*(Mi-x));
    z2=max((sum(wi.*x)/W)-1,0);
    %z2=sum(wi.*x);
    z=[z1,z2]';
    
    W1=sum(wi.*x);
    V1=sum(vi.*x);
    
    sol.V1=V1;
    sol.W1=W1;
    vol=1;
end