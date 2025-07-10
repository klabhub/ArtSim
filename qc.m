function [r,mse] = qc(x,truth,tacsFreq,sf)
%function [r,mse] = qc(x,truth)
% Helper function to quantify the match between the cleaned signal and ground truth
% Using Pearson correlation and mean-squared error.
% INPUT
% x -  Cleaned signal
% truth - Ground truth signal
% tacsFreq - The frequency that will be removed from the x and truth before
%               comparison [NaN]
% sf - Sampling frequency of x and truth  [NaN].
% OUTPUT
% r = Pearson correlation between signal and truth
% mse = Mean squared error difference between signal and truth  in muV.

arguments
    x double 
    truth double
    tacsFreq (1,1) double = NaN
    sf (1,1) double = NaN
end


x= x(:);
truth=truth(:);
if ~isnan(tacsFreq)
    % Remove the tacsFreq from consideration
    x= FBAR(x,tacsFreq,sf);
    truth = FBAR(truth,tacsFreq,sf);    
end
mse = sqrt(median((x-truth).^2,'omitnan'))/1e-6;
out= isnan(x) | isnan(truth);
r = corr(x(~out),truth(~out),'type','Pearson');
end

