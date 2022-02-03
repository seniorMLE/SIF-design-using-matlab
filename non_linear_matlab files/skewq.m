function [ A ] = skewq( B )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

A=[   0   -B(1) -B(2) -B(3);
     B(1)   0    B(3) -B(2);
     B(2) -B(3)   0    B(1);
     B(3)  B(2) -B(1)   0];
end
