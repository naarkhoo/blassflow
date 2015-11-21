function [ probability_of_x ] = expo_likelihood( lambda, x )
% calculate the probability of x, random number form exp(lambda)

    probability_of_x = lambda * exp(- lambda * x);

end

