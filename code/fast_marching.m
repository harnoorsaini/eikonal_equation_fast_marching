function travelTime = fast_marching(sources,domainBounds,plot,delta,dim, ...
    domainLength)

% -------------------------------------------------------------------------
% XXX
% types
% 0 - alive
% 1 - narrow band
% 2 - far away
%delta=1;

numSources = size(sources,1);

% -------------------------------------------------------------------------
% Define neighbour connectivity; currently: 
%
%         [3](0,1)
%              |
%              |
% [2](-1,0)--(0,0)--[1](1,0)
%              | 
%              |
%         [4](0,-1)
%

neighbourStencil = [1 0;-1 0;0 1;0 -1];
numNeighbours = size(neighbourStencil,1);

% -------------------------------------------------------------------------
% set up empty values at FD grid points
% speed values mm/s
nodeSpeed = ones(domainBounds(1),domainBounds(2));
%nodeSpeed(1:4,1:4) = zeros;
%nodeSpeed(4:6,3) = zeros;
% travel times - to be solved for; initiall inifinity 
travelTime = ones(domainBounds(1),domainBounds(2))*inf;
% set state: frozen(0); narrow band (1); or unknown (2)
nodeState = 2*ones(domainBounds(1),domainBounds(2));

% set source points to frozen
% sub2ind to convert 2d index location to 1d vector location
idx_Source = sub2ind(domainBounds,sources(:,1),sources(:,2)); 
travelTime(idx_Source) = 0;
nodeState(idx_Source) = 0;

neighbourNodes(numNeighbours,dim) = zeros;

% -------------------------------------------------------------------------
% define narrow band for the number of sources by finding the neighbours
% around all the sources
for idx_source = 1:numSources
	for idx_dim = 1:dim
        % use stencil at each source
        neighbourNodes(:,idx_dim) = ... 
            sources(idx_source,idx_dim)+neighbourStencil(:,idx_dim);
        
        % check if any neighbours lie outside the FD grid
        Logidx =  find( neighbourNodes(:,idx_dim)>domainBounds(idx_dim) ... 
            | neighbourNodes(:,idx_dim)<0) ;
        
        % use logical indexing to ommit these "ghost neighbours" 
        neighbourNodes(Logidx,idx_dim) = sources(idx_source,idx_dim);
    end
    % corresponding linear indices
    idx_Neighbour = ...
        sub2ind(domainBounds,neighbourNodes(:,1),neighbourNodes(:,2));
    
    nodeState(idx_Neighbour) = (nodeState(idx_Neighbour)~=0);
    
    % find the travel time using first order backward differences 
    travelTime(idx_Neighbour)=delta./nodeSpeed(idx_Neighbour);
end

% -------------------------------------------------------------------------
stop=0;
iter=1;

% local travel time computation
travelTime_unkown(dim,1) = zeros;

% local travel times storage
travelTime_neighb(dim,1) = zeros;

coordNeighb_Neighb(numNeighbours,dim) = zeros;

while ~stop,

	% continue till all nodes have been traversed
	lin_narrowBand = find(nodeState==1); 
    % ---------------------------------------------------------------------
	if ~isempty(lin_narrowBand),
        % find the coordinate of the smallest T
		[~,idx_neighMinTravelTime] = min(travelTime(lin_narrowBand));
		
        idx_neighMinT = lin_narrowBand(idx_neighMinTravelTime);
		
        % freeze this point
		nodeState(idx_neighMinT) = 0;
        
		% find FD grid location for neighbour with lowest traven time
		[nodeCoord_x,nodeCoord_y]=ind2sub(domainBounds,idx_neighMinT);
		
        currNode=[nodeCoord_x,nodeCoord_y];

        % find the coordinates of the neighbours to the new frozen point
		for idx_dim=1:dim
			neighbourNodes(:,idx_dim) = ...
                currNode(idx_dim)+neighbourStencil(:,idx_dim);
			
            Logidx =find(neighbourNodes(:,idx_dim)>domainBounds(idx_dim) ...
                | neighbourNodes(:,idx_dim)<=0);
            			
			neighbourNodes(Logidx,idx_dim) = currNode(idx_dim) ... 
                -(neighbourNodes(Logidx,idx_dim)-currNode(idx_dim));
        end

        % isolate only those neighbours which are unknown 
        
        % find linear index of all new neighbours
		lin_idx_Neighbours = sub2ind(domainBounds,...
            neighbourNodes(:,1),neighbourNodes(:,2));

        % identify index of those neighbours which are unknown
        Logidx_unknownNeighb = nodeState(lin_idx_Neighbours)==2;
        
        % find the FD nodal index of these unknown neighbours
        lin_idx_NewNeighb = lin_idx_Neighbours(Logidx_unknownNeighb);
		      
        % invite the unknown neighbours!
        nodeState(lin_idx_NewNeighb) = 1;
        
		for idx_Neighb=1:numNeighbours
			currentNeighbour = neighbourNodes(idx_Neighb,:);

            NodeSpeed_curr = nodeSpeed(currentNeighbour(1), ... 
                currentNeighbour(2));

            % find the coordinates of the neighbours neighbours
			for idx_dim=1:dim
				coordNeighb_Neighb(:,idx_dim) = ... 
                    currentNeighbour(idx_dim) + neighbourStencil(:,idx_dim);
				
                Logidx = find(coordNeighb_Neighb(:,idx_dim) ... 
                    >domainBounds(idx_dim) | ... 
                    coordNeighb_Neighb(:,idx_dim)<=0);
                
				coordNeighb_Neighb(Logidx,idx_dim) =  ... 
                    currentNeighbour(idx_dim) - ... 
                    (coordNeighb_Neighb(Logidx,idx_dim) ... 
                    - currentNeighbour(idx_dim));
			end

            % find the smallest travel time of the neighbours neighbour
            % bascially the values will either be infinity or some finite
            % value on either side of the neighour
			for idx_dim=1:dim
                % minimum in x-direction
				temp_1 = travelTime(coordNeighb_Neighb(2*idx_dim-1,1)...
                    ,coordNeighb_Neighb(2*idx_dim-1,2));
				
                % minimum in y-direction
                temp_2 = travelTime(coordNeighb_Neighb(2*idx_dim,1)...
                    ,coordNeighb_Neighb(2*idx_dim,2));
				
                travelTime_neighb(idx_dim) = min(temp_1,temp_2);
            end

            % sort known travel times according to lowest travel times
			[travelTime_neighb,~] = sort(travelTime_neighb);
            
            % in nD, the nth sorted travel time is the slowest
			idx_slowerT = dim; 
            
            % find the minimum travel time that can occur (to a "linear"
            % neighbour
			travelTime_unkown(1) = travelTime_neighb(1) ... 
                + delta/NodeSpeed_curr;
			travelTime_unkown(idx_slowerT) = travelTime_unkown(1);
			cont = 0;
            
			% while the fastest travel time is slower the slowest travel
			% time of the neighbour 
			while travelTime_unkown(idx_slowerT-1) ...
                    > travelTime_neighb(idx_slowerT) && ~cont 

                % solve the quadratic equation for the "corner" neigbhours
                alpha = 2;
                beta = -2*sum(travelTime_neighb(1:idx_slowerT));
                gamma = travelTime_neighb(1)^2+travelTime_neighb(2)^2 ...
                    - delta^2/NodeSpeed_curr^2;
                T_tmp1 = (-beta + sqrt(beta^2-4*alpha*gamma) ) / (2*alpha);
                T_tmp2 = (-beta - sqrt(beta^2-4*alpha*gamma) ) / (2*alpha);
                
                % take the slower travel time
                travelTime_unkown(idx_slowerT) = max(T_tmp1,T_tmp2);
				
				if idx_slowerT < dim-2
					idx_slowerT = idx_slowerT + 1;
				else
					cont = 1;
				end
			end

			%...
			if nodeState(currentNeighbour(1),currentNeighbour(2)) == 1
				travelTime(currentNeighbour(1),currentNeighbour(2)) = travelTime_unkown(idx_slowerT);
            end
		end

	else
		stop = 1;
    end
    
	iter = iter + 1;
    
	%fprintf('paso: %d \n',iter);
	if (plot)
		figure(1);
		mesh([0:delta:domainLength(1)],[0:delta:domainLength(2)],travelTime);
		%axis([0 domainLength(1)+1 0 domainLength(2)+1 0 domainLength(2)+1])
		drawnow;
	end
	%colormap(gray(256));
	%     figure(2);
	%     mapa_clas = [1 1 1;hsv(3)];
	%     imshow(tipos+1,mapa_clas);
	%    h=gcf;
	%    set(h,'Position',[148 104 528 528]);
	%    hh=get(h,'Children')
	%    set(hh,'Position',[0.1 0.1 0.8 0.8]);
	%    pause;
end

%if (plot)
%	figure
	%image(255./(max(max(T))-min(min(T)))*(T-min(min(T))));
	%colormap(gray(256))
%	mesh(travelTime);
%end