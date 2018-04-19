function X = CreateGriddedGaussian( D, ptsperdim, center, var, pnorm );

if nargin<5 || isempty(pnorm),    pnorm = 2;  end

% Create grid in ptsperdim^D
str_eval = '[';
for k = 1:D,
    str_eval = [str_eval,sprintf('X_%d,',k)];
end;
str_eval(end) = ']';
str_eval = [str_eval,'=ndgrid(ptsperdim);'];
eval(str_eval);

% Construct Gaussians
str_eval    = 'exp((';
for k = 1:D,
    str_eval    = [str_eval sprintf('-abs(X_%d-center(%d)).^%d',k,k,pnorm)];
end;
eval(sprintf('X = %s )/(%f^(%d/2)));',str_eval, var,pnorm));

return;
