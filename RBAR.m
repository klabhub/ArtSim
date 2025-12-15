function v = RBAR(v,tacsArtifact)
% Use a linear regression to remove the artifact.
% Solve this regression for each segment
% vRecorded = offset + slope*tacsArtifact
% then remove this from the recorded signal by substraction 
% vClean = vRecorded - slope*tacsArtifact
% Note that the offset is kept to avoid removng the mean in each
% segment, which could remove signal for short segments.
arguments
    v (:,:) double    %[nrSamplesPerSegment nrSegments]
    tacsArtifact (:,:) double %[nrSamplesPerSegment nrSegments]
end

[nrSamplesPerSegment,nrSegments] = size(v); %#ok<ASGLU>
beta  = zeros(2,nrSegments);
for i=1:size(v,2)
    stay = ~isnan(v(:,i)) & ~isnan(tacsArtifact(:,i));  % Keep non-nan
    D = [ones(sum(stay),1) tacsArtifact(stay,i)];
    beta(:,i) = D \ v(stay, i); % Regress    
end
v = v - beta(2,:).*tacsArtifact; % Remove signal that matches the tACS artifact
end
