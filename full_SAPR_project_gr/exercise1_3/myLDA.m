function A = myLDA(Samples, Labels, NewDim)
% Input:    
%   Samples: The Data Samples 
%   Labels: The labels that correspond to the Samples
%   NewDim: The New Dimension of the Feature Vector after applying LDA
    
	[NumSamples NumFeatures] = size(Samples);
    
    A=zeros(NumFeatures,NewDim);

    NumLabels = length(Labels);
    if(NumSamples ~= NumLabels) then
        fprintf('\nNumber of Samples are not the same with the Number of Labels.\n\n');
        exit
    end
    Classes = unique(Labels);
    NumClasses = length(Classes);  %The number of classes
    
    P = zeros(1,NumClasses);
    Sw = zeros(NumFeatures);
    Sb = zeros(NumFeatures); 

    %predefine variable to hold class mean
    mu = zeros(NumClasses , NumFeatures);

    %For each class i
	%Find the necessary statistics
    for i=1:NumClasses 
        %Calculate the Class Prior Probability
	    P(i)=sum(Labels==(i-1))/NumSamples;
      
        %Calculate the Class Mean 
	    mu(i,:)=mean(Samples( Labels == (i-1), :));
        
        %Calculate the Global Mean
	    m0=mean(Samples);

        %Calculate the Within Class Scatter Matrix
	    Sw= Sw + P(i)*cov(Samples(Labels==(i-1),:));

        %Calculate the Between Class Scatter Matrix
	    Sb=Sb + P(i)*(mu(i,:) - m0)'*(mu(i,:) - m0);  

        %Eigen matrix EigMat=inv(Sw)*Sb
        EigMat = inv(Sw)*Sb;

    end 
   
    %Perform Eigendecomposition
    [V,D]=eig(EigMat);
    eigenval= diag(D); % Eigenvalues lie in the diagonal of the covariance matrix
    eigenvec = V; %Corresponding eigenvectors
    [eigenval, order] = sort(eigenval, 'descend');  %Sorting eigenvalues and returning the order
    eigenvec = eigenvec(:, order);
    
    %Select the NewDim eigenvectors corresponding to the top NewDim
    %eigenvalues (Assuming they are NewDim<=NumClasses-1)
	%% You need to return the following variable correctly.
	A = eigenvec(:, 1:NewDim);  % Return the LDA projection vectors
