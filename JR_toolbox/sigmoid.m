function y = sigmoid(x,phase,slope)
% y = sigmoid(x,phase,slope)

y=1./(1+phase*exp(-slope*(x)));
