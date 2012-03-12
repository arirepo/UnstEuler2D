%this script tests the compact schemes for integration which we excluded
%the initial function value f1 by slicing the matrices A and B
clear all
clc

%format long

n = 5;
[A B f1] = S(n,2,[3 3],2);
%A
%B
A\f1
SS = A\B
