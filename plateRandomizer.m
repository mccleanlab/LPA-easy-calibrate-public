clear all; close all; clc;

%% Generates random map of wells to measure

numWells = 8;
x = randperm(24);
x = x(1:numWells);
plate = [(1:6); (7:12); (13:18); (19:24)];
[v, i] = intersect(plate,x);
plate(i) = -plate(i);
plate(plate>0)=0; 
plate(plate<0)=1;
disp(plate)
