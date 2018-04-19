function K = GenerateIntegralOperatorFromSets( X,Y,Opts )

% Uses: IPDM, may be downloaded at 
% http://www.mathworks.com/matlabcentral/fileexchange/18937

distances = ipdm(X,Y);

K = 1./distances;

return;