function y=sumOfTwoGaussians(beta,x)

beta=max(beta,[0.1 2 0.5 0.1 7 0.5]);
y=beta(1)*normpdf(x,beta(2),beta(3)) + beta(4)*normpdf(x,beta(5),beta(6));
y=y(:);