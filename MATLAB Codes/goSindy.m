function outp = goSindy(depVar,Theta,lambdaGuess,nlambda,lambda,plt)

if isempty(lambda)
    lambdaRefN0 = Inf;
    lambdaRef0 = -Inf;
    ratioN0 = lambdaRefN0/lambdaGuess;
    ratio0 = lambdaRef0/lambdaGuess;
    
    while ratioN0 > 1.5 || ratio0 < 0.5
        Xi = sparsifyDynamics(Theta,depVar,lambdaGuess);
        NZT = sum(abs(Xi)>eps);
        if NZT ~= 0
            ratioN0 = lambdaRefN0/lambdaGuess;
            lambdaRefN0 = lambdaGuess;
            lambdaGuess = lambdaGuess*1.5;
        elseif NZT == 0
            ratio0 = lambdaRef0/lambdaGuess;
            lambdaRef0 = lambdaGuess;
            lambdaGuess = lambdaGuess*0.5;
        end
    end
    
    lambdaMax = max([lambdaRef0,lambdaRefN0,lambdaGuess]);
    lambda = linspace(0,lambdaMax,nlambda);
else
    if length(lambda)>1
        error("lambda must be a scalar input.")
    end
    nlambda = 1;
end

Xi = zeros(size(Theta,2),nlambda);
idxNonZero = zeros(size(Theta,2),nlambda);
MSE = zeros(nlambda,1);
AIC = zeros(nlambda,1);
depVarEst = zeros(length(depVar),nlambda);

for i = 1:length(lambda)
    Xi(:,i) = sparsifyDynamics(Theta,depVar,lambda(i));
    depVarEst(:,i) = Theta*Xi(:,i);
    idxNonZero(:,i) = (abs(Xi(:,i))>eps);
    MSE(i) = (1/length(depVar))*(norm(depVar-depVarEst(:,i))^2);
    AIC(i) = length(depVar)*log(MSE(i))+2*sum(idxNonZero(:,i))...
        +(2*sum(idxNonZero(:,i))*(sum(idxNonZero(:,i))+1))...
        /(length(depVar)-sum(idxNonZero(:,i))-1);
end

if strcmpi(plt,'y')
    figure(100)
    subplot(3,1,1)    
    plot(lambda,MSE,'.-','LineWidth',2,'MarkerSize',15)
    xlabel('\lambda')
    ylabel('MSE')
    set(gca,'LineWidth',1.5,'FontSize',12)
    subplot(3,1,2)
    plot(lambda,AIC,'.-','LineWidth',2,'MarkerSize',15)
    xlabel('\lambda')
    ylabel('AIC')   
    set(gca,'LineWidth',1.5,'FontSize',12)
    subplot(3,1,3)
    plot(lambda,sum(idxNonZero,1),'.-','LineWidth',2,'MarkerSize',15)
    xlabel('\lambda')
    ylabel('# of terms')
    set(gca,'LineWidth',1.5,'FontSize',12)
end

outp.lambda = lambda;
outp.MSE = MSE;
outp.Theta = Theta;
outp.Xi = Xi;
outp.idxNonzero = idxNonZero;
outp.NZT = sum(idxNonZero,1);
outp.depVarEst = depVarEst;

end

%%-------------------------------------------------------------------------

function Xi = sparsifyDynamics(Theta,dXdt,lambda)

% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% compute Sparse regression: sequential least squares
Xi = Theta\dXdt;  % initial guess: Least-squares

%******
ndXdt=size(dXdt,2);

% lambda is our sparsification knob.
for k=1:10
    idxSmall = (abs(Xi)<lambda);   % find small coefficients
    Xi(idxSmall)=0;                % and threshold
    for ind = 1:ndXdt                   % n is state dimension
        idxBig = ~idxSmall(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(idxBig,ind) = Theta(:,idxBig)\dXdt(:,ind); 
    end
end

end