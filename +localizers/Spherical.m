function [ y ] = Sperical( rho, theta1, theta2 )
%SPERICAL Summary of this function goes here
%   Detailed explanation goes here

y = arrayfun(@(rhoi) scalarSperical(rhoi,theta1, theta2), rho);
end


function y = scalarSperical(rho, theta1, theta2)
if (abs(rho) < theta2)
    y = theta1 * (1 - 2/pi * (abs(rho)/theta2 * sqrt(1 - (abs(rho)/theta2)^2) + sin(abs(rho)/theta2)^(-1)));
else
    y = 0;
end
end