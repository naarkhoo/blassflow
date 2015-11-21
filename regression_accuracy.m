function [ err ] = error( z1, z, method )
%RMS_ERROR Summary of this function goes here
%   Detailed explanation goes here

switch(method)
    case 'mse'
        err = mean(mean((z1 - z).^2));
    case 'rmse'
        err = sqrt(mean((z1 - z).^2));
end
        
end

