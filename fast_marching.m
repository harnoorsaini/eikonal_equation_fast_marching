function travelTime = fast_marching(sources,domainBounds,plot,delta,dim)

% -------------------------------------------------------------------------
% XXX
% types
% 0 - alive
% 1 - narrow band
% 2 - far away
%delta=1;

numSources = size(sources,1);
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

% set up empty values at FD grid points
% speed values mm/s
nodeSpeed = ones(domainBounds(1),domainBounds(2));
nodeSpeed(14,5:size(nodeSpeed,1)-5) = 0;
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

% define narrow band for the number of sources by finding the neighbours
% around all the sources
for idx_source = 1:numSources
	for idx_dim = 1:dim
        % use stencil at each source
        neighbourNodes(:,idx_dim) = ... 
            sources(idx_source,idx_dim)+neighbourStencil(:,idx_dim);
        
        % check if any neighbours lie outside the FD grid
        Logidx =  neighbourNodes(:,idx_dim)>domainBounds(idx_dim) ... 
            | neighbourNodes(:,idx_dim)<0 ;
        
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


stop=0;
iter=1;

% local travel time computation
travelTime_local(dim,1) = zeros;

% local travel times storage
m(dim,1) = zeros;

coordNeighb_Neighb(numNeighbours,dim) = zeros;

while ~stop,

	% find neighbour with lowest travel time in the narrow band.
	lin_narrowBand = find(nodeState==1); % returns a linear index
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
			
            Logidx= neighbourNodes(:,idx_dim)>domainBounds(idx_dim) ...
                | neighbourNodes(:,idx_dim)<=0 ;
            			
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
				
                Logidx = coordNeighb_Neighb(:,idx_dim) ... 
                    >domainBounds(idx_dim) | ... 
                    coordNeighb_Neighb(:,idx_dim)<=0;
                
				coordNeighb_Neighb(Logidx,idx_dim) =  ... 
                    currentNeighbour(idx_dim) - ... 
                    (coordNeighb_Neighb(Logidx,idx_dim) ... 
                    - currentNeighbour(idx_dim));
			end

            % find the smallest travel time of the neighbours neighbour
			for idx_dim=1:dim
                % minimum in x-direction
				temp_1 = travelTime(coordNeighb_Neighb(2*idx_dim-1,1)...
                    ,coordNeighb_Neighb(2*idx_dim-1,2));
				
                % minimum in y-direction
                temp_2=travelTime(coordNeighb_Neighb(2*idx_dim,1)...
                    ,coordNeighb_Neighb(2*idx_dim,2));
				
                m(idx_dim)=min(temp_1,temp_2);
			end

            
            % for a corner node; m(1) and m(2) will both have non-infinite
            % values; but for a "linear" node; either m(1) or m(2) will be
            % infinite (as the wave has not reached the linear neighbours
            % yet) - WELL FOR ***ONE*** SOURCE AT LEAST!

            % sort travel times according to lowest travel times
			[m,ind_m]=sort(m);
            % this is hard coded? (dimensionality?)
			idx_source=2; % always the "slower" arrival time 
            % travel time of neigbhour + step-size/speed; this is the
            % minimal travel time there can be (note using m(1); which has
            % already been sorted) 
			travelTime_local(1)=m(1)+delta/NodeSpeed_curr;
			travelTime_local(idx_source)=travelTime_local(1);
			me_pase_de_dim=0;
			% we have to check this account, it gives sweet potatoes
            % maybe a check for being next to another source?
			while travelTime_local(idx_source-1)>m(idx_source) & ~me_pase_de_dim
				%fprintf('8');
                % solving the roots of the quadratic equation for corner
                % points
				%u(k)=sum(m(1:k))+sqrt( sum(m(1:k))^2 + k*delta^2/F_curr + k*sum(m(1:k).^2) );
				%u(k)=sum(m(1:k))+sqrt( sum(m(1:k))^2 + k*delta^2/F_curr - k*sum(m(1:k).^2) );
                alpha = 2;
                beta = -2*sum(m(1:idx_source));
                gamma = m(1)^2+m(2)^2 - delta^2/NodeSpeed_curr^2;
                u_tmp1 = (-beta + sqrt(beta^2-4*alpha*gamma) ) / (2*alpha);
                u_tmp2 = (-beta + sqrt(beta^2-4*alpha*gamma) ) / (2*alpha);
                travelTime_local(idx_source) = max(u_tmp1,u_tmp2);
				%u(k)=u(k)/k;
				if idx_source<dim-2
					idx_source=idx_source+1;
				else
					me_pase_de_dim=1;
				end
			end

			%T(current_vecino(1),current_vecino(2))=u(k);
			% REVIEW: this may be wrong
			% I'm changing the value even if it's not NW I'd have to check this.
			% able that there is another better place to do it...
			% k or k-1??????
			if nodeState(currentNeighbour(1),currentNeighbour(2))==1
				travelTime(currentNeighbour(1),currentNeighbour(2))=travelTime_local(idx_source);
			end
			%fprintf('9');
			% keyboard
		end

	else
		stop=1;
	end
	iter=iter+1;
	%fprintf('paso: %d \n',iter);
	if (plot)
		figure(1);
		mesh(travelTime);
		axis([0 50 0 50 0 50])
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

if (plot)
	figure
	%image(255./(max(max(T))-min(min(T)))*(T-min(min(T))));
	%colormap(gray(256))
	mesh(travelTime);
end