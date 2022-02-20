function model=CreateModel()
    
    wi=[19 10 10 22 40 24 26 27 16 12];
    vi=[66 70 39 45 41 28 71 20 63 73];
    Mi=[14 14 14 20 16 19 10 14 13 14];
    W=1548;
    
    model.W=W;
    model.wi=wi;
    model.vi=vi;
    model.Mi=Mi;
    model.n=numel(wi);
end