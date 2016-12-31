function geoDist = GeodesicSaliency(adjcMatrix, bgIds, colDistM, clip_value,bg_pro)

spNum = size(adjcMatrix, 1);

% Calculate pair-wise geodesic distance
adjcMatrix_lb = adjcMatrix;%adjacent matrix with boundary SPs linked
adjcMatrix_lb(bgIds, bgIds) = 1;
 
[row,col] = find(adjcMatrix_lb);

% Here we add a virtual background node which is linked to all background
% super-pixels with 0-cost. To do this, we padding an extra row and column 
% to adjcMatrix_lb, and get adjcMatrix_virtual.
adjcMatrix_virtual = full(sparse([row; repmat(spNum + 1, [length(bgIds), 1]); bgIds], ...
    [col; bgIds; repmat(spNum + 1, [length(bgIds), 1])], 1, spNum + 1, spNum + 1));
% Specify edge weights for the new graph
colDistM_virtual = zeros(spNum+1);
colDistM_virtual(1:spNum, 1:spNum) = colDistM;
adjcMatrix_virtual = tril(adjcMatrix_virtual, -1);

colDistM_virtual = zeros(spNum+1);
colDistM_virtual(1:spNum, 1:spNum) = colDistM;
colDistM_virtual = max(0, colDistM_virtual - clip_value);
colDistM_virtual = colDistM_virtual/max(colDistM_virtual(:));
colDistM_virtual(end,bgIds) = bg_pro(:);
colDistM_virtual(bgIds,end) = bg_pro(:);
% Specify edge weights for the new graph


edgeWeight = colDistM_virtual(adjcMatrix_virtual > 0);
% edgeWeight = max(0, edgeWeight - clip_value);
% edgeWeight = edgeWeight/max(edgeWeight(:));
% edgeWeight(end,bgIds) = bg_pro(:);
% edgeWeight(bgIds,end) = bg_pro(:);
geoDist = graphshortestpath(sparse(adjcMatrix_virtual), spNum + 1, 'directed', false, 'Weights', edgeWeight);
geoDist = geoDist(1:end-1); % exclude the virtual background node
%geoDist = geoDist / max(geoDist(:)); 


