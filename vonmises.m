function v = vonmises(x,theta,kappa)
% function v = vonmises(x,theta,kappa)
% Returns the von Mises function at x for a given preferred
% direction (theta) and concentration parameter (kappa)
arguments
    x (:,:) double
    theta (1,1) double
    kappa  (1,1) double
end

v= (exp(kappa*cos(x-theta)) ./ (2*pi*besseli(0,kappa)));

