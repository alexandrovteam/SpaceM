function vSignal = UnPatchifySignal( cPatches, cWindowSize, cTau )

%
% Splits cSignal into pieces
%

if nargin<3,
    cTau = 1;
end;

lN = size(cPatches,1);

vSignal = zeros(lN,1);

for k = 1:lN,
    vSignal(k) = cPatches(k,floor(cWindowSize/2));
end;

return;