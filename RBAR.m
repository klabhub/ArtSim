function v = RBAR(v,tacsSignal)
% Use a linear regression to remove the artifact.
% Solve this regression
% vRecorded = slope*tacsSignal
% then remove 
% vClean = vRecorded-slope*tacsSignal
arguments
    v (:,:) double
    tacsSignal (:,:) double
end

slope  = zeros(1,size(v,2));
for i=1:size(v,2)
    stay = ~isnan(v(:,i)) & ~isnan(tacsSignal(:,i));
    slope(i) = tacsSignal(stay,i)\v(stay,i);
end
v = v-slope.*tacsSignal;
end
