function v = RBAR(v,tacsArtifact)
% Use a linear regression to remove the artifact.
% Solve this regression for each segment
% vRecorded = offset + slope*tacsArtifact
% then remove this from the recorded signal by substraction 
% vClean = vRecorded-offset - slope*tacsArtifact
arguments
    v (:,:) double    %[nrSamplesPerSegment nrSegments]
    tacsArtifact (:,:) double %[nrSamplesPerSegment nrSegments]
end

[nrSamplesPerSegment,nrSegments] = size(v); %#ok<ASGLU>
beta  = zeros(2,nrSegments);
for i=1:size(v,2)
    stay = ~isnan(v(:,i)) & ~isnan(tacsArtifact(:,i));  % Keep non-nan
    D = [ones(sum(stay),1) tacsArtifact(stay,i)];
    beta(:,i) = regress(v(stay,i),D(stay,:)); % Regress    
end
v = v - 0*beta(1,:) - beta(2,:).*tacsArtifact; % Remove signal that matches the tACS artifact (plus the mean of the segment)
end
