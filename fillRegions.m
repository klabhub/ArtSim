function [filled,nrInRegion,startRegionIx,stopRegionIx] = fillRegions(vec, fillIn)
% For a vector of logicals, fill in regions with the value true if the
% number of elements between them is less than fillIn (in Hz). 
% Returns the filled vector, the number of elements in each region and the
% start and stop indices of all the true-regions.
arguments
    vec (:,:) logical 
    fillIn (1,1) double {mustBeNonnegative} 
end

% Initialize the output vector as the input vector
vec = logical(vec(:));
filled = vec;
% Find the starting indices of regions with 1's in the vector
startRegionIx = 1+find(diff(vec)>0);
stopRegionIx =  find(diff(vec)<0);
if vec(1)
    startRegionIx = [1;startRegionIx];
end
if vec(end)
    stopRegionIx = [stopRegionIx;numel(vec)];
end
for regionCntr = 1:length(startRegionIx)-1
    if (startRegionIx(regionCntr+1)-stopRegionIx(regionCntr)-1)<= fillIn
        filled((stopRegionIx(regionCntr)+1):(startRegionIx(regionCntr+1)-1))=true;
    end
end
startRegionIx = 1+find(diff(filled)>0);
stopRegionIx =  find(diff(filled)<0);
if filled(1)
    startRegionIx = [1;startRegionIx];
end
if filled(end)
    stopRegionIx = [stopRegionIx;numel(vec)];
end
nrInRegion  = stopRegionIx- startRegionIx+1;

end
