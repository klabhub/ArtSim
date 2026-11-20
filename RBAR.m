function v = RBAR(v,tacsArtifact)
% Use a linear regression to remove the artifact.
% Solve this regression for each segment
% vRecorded = offset + slope*(tacsArtifact-artifactBaseline)
% then remove the time-varying artifact from the recorded signal.
% The baseline is estimated once across the complete proxy. This encodes
% the assumption that constant current produces no artifact, while changes
% in the measured stimulation voltage represent artifact.
arguments
    v (:,:) double    %[nrSamplesPerSegment nrSegments]
    tacsArtifact (:,:) double %[nrSamplesPerSegment nrSegments]
end

[nrSamplesPerSegment,nrSegments] = size(v); %#ok<ASGLU>

validArtifact = ~isnan(tacsArtifact) & ~isinf(tacsArtifact);
artifactBaseline = median(tacsArtifact(validArtifact));
if isempty(artifactBaseline) || ~isfinite(artifactBaseline)
    error('RBAR:InvalidArtifactProxy', ...
        'The artifact proxy contains no finite samples.');
end
% Constant current/voltage should produce no artifacts; remove 
artifactDeviation = tacsArtifact - artifactBaseline;

beta  = zeros(2,nrSegments);
for i=1:size(v,2)
    stay = ~isnan(v(:,i)) & ~isinf(v(:,i)) & validArtifact(:,i);
    D = [ones(sum(stay),1) artifactDeviation(stay,i)];
    beta(:,i) = D \ v(stay, i); % Regress    
end
v = v - beta(2,:).*artifactDeviation; % Remove time-varying artifact
end
