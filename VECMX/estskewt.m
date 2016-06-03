function [x,fval,exitflag]=estskewt(data)

%% Objective function
    function [LL, likelihoods] = skewt_llf(pars, data)
nu = pars(1); %Degrees of Freedom parameter
lambda = pars(2); %Skewness parameter
[T,k] = size(data);

logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
c = exp(logc);
a = 4*lambda.*c.*((nu-2)./(nu-1));
logb = 0.5*log(1 + 3*lambda.^2 - a.^2);
b = exp(logb);

find1 = (data<(-a./b));
find2 = (data>=(-a./b));
LL1 = logb + logc - (nu+1)/2.*log(1+1./(nu-2).*((b.*data(find1)+a)./(1-lambda)).^2);
LL2 = logb + logc - (nu+1)/2.*log(1+1./(nu-2).*((b.*data(find2)+a)./(1+lambda)).^2);
LL = sum(LL1) + sum(LL2);
LL = -LL;

if nargout>1
     likelihoods=zeros(size(data));
     likelihoods(find1)=LL1;
     likelihoods(find2)=LL2;
end
    end
%% Optimisation Parameters & Starting Values
nu=5;
lambda=0;
pars=[nu;lambda];
parsL=[0;-15];
parsU=[Inf;15];
% Initial estimates
x0 = [5;0];
%% Optimisation Options Set - fmincon
optopt=optimset('Algorithm', 'interior-point', 'MaxIter',100,'Display','Iter','LargeScale','on')

%% Optimisation Subroutine
[x,fval,exitflag]=fmincon(@skewt_llf,pars, [],[],[],[],parsL,parsU,[],optopt,data);
pars=x;
[LL,likelihoods]=skewt_llf(pars,data);

end