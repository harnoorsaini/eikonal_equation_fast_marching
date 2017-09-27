% =========================================================================
% 2D Fast Marching Algorithm to solve the Eikonal Equation
% =========================================================================
% Harnoor Saini 
% Sep-2017
%
% Continuum Biomechanics and Mechanobiology 
% University of Stuttgart
% Stuttgart, Germany
%
% -------------------------------------------------------------------------
% Original script "2D fast marching algorithm" by Juan Cardelino (2013) can
% be found at https://de.mathworks.com/matlabcentral/fileexchange/...
% ...18529-2d-fast-marching-algorithm
% -------------------------------------------------------------------------
% Main changes:
%
% 1. Core finite difference approximation (for corner nodes) has
%    been modified from the original version
%
% 2. Commeting, formatting (conversion to english variables) and general
%    clean up of code
%
% 3. Generalisation of code and function I/O
%
% 4. Speed up
%
% =========================================================================
% TO DO:
% - make find_neighbours into a seperate function
% - a way to nicely/quickly define the speeds (including obstacles)
%
%
%
%
% -------------------------------------------------------------------------
% RESET WORKSPACE
clear 
close all
clc

% -------------------------------------------------------------------------
% USER INPUTS

% dimensionality (only works for 2!)
dim = 2;

% plotting
plotting = 1;

% set finite difference grid geometry [x y]
domainLength = [10 10]; %mm
delta = 1; %mm

% set source locations in mm [x y] @ 1 per row
randSources = 0;
if randSources
    sources = [domainLength(1)*rand domainLength(2)*rand; ... 
    domainLength(1)*rand domainLength(2)*rand; ...
    domainLength(1)*rand domainLength(2)*rand; ...
    domainLength(1)*rand domainLength(2)*rand]; %mm
else
    sources = [5 5];
end
numSources = size(sources,1);

% -------------------------------------------------------------------------
% SET UP FINITE DIFFERENCE GRID

idx_x = 1;
idx_y = 2;
nodeCoords_x = 1:delta:domainLength(idx_x); 
nodeCoords_y = 1:delta:domainLength(idx_y);

numNodes_x = size(nodeCoords_x,2)+1;
numNodes_y = size(nodeCoords_y,2)+1;
domainNodesBounds = [numNodes_x,numNodes_y];
% sources
sourcesNodes(numSources,idx_y) = zeros;
for idx_source = 1:numSources
    sourcesNodes(idx_source,:) = round(domainNodesBounds.*...
        sources(idx_source,:)./domainLength);
end

% -------------------------------------------------------------------------
% PERFORM MAST MARCHING TO FIND TRAVEL TIMES

T = fast_marching(sourcesNodes,domainNodesBounds,plotting,delta,dim, ...
    domainLength);

%if plotting
%    mesh(T);
%end



