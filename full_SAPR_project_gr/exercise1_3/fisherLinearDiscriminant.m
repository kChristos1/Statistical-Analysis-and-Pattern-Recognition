function v = fisherLinearDiscriminant(X1, X2)

    [m,n] = size(X1);

    m1 = size(X1, 1);
    m2 = size(X2, 1);

    mu1 = mean(X1); % mean value of X1
    mu2 = mean(X2); % mean value of X2


    X1_c = X1 - mu1;
    X2_c = X2 - mu2; 
 

    S1 = (X1_c.' * X1_c);  % scatter matrix of X1
    S2 = (X2_c.' * X2_c); % scatter matrix of X2


    Sw = (S1+S2)/2; % Within class scatter matrix  
    invS=inv(Sw);
   
    v = invS*(mu1-mu2).';
   
    
    v = v/(norm(v)); % return a vector of unit norm
