% SUMMARY:
% This script randomly generates a map showing the subset of wells of a
% 24-well plate that should be measured when determining the relationship
% between IRIS value and light intensity output.
%  
% INPUTS: 
% None
%  
% OUTPUTS:
% plateMap: binary map showing wells to be measured (1) or ignored (0)
%
% WRITTEN BY:
%   Kieran Sweeney
%   University of Wisconsin-Madison
%   Department of Biomedical Engineering
%   1550 Engineering Drive ECB 3156
%   Madison, WI 53705
%   ksweeney2@wisc.edu
%
% Last revised on August 30, 2018

%% Prepare to run script
clearvars; close all; clc;

%% Set number of wells to be measured
numWells = 8;

%% Generate random map of wells to be measured
x = randperm(24,numWells);
plateMap = [(1:6); (7:12); (13:18); (19:24)];
[v, i] = intersect(plateMap,x);% 
plateMap(i) = -plateMap(i);
plateMap(plateMap>0)=0; 
plateMap(plateMap<0)=1;
disp(plateMap)
