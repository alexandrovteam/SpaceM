function vPatches = PatchifySignal( cSignal, cWindowSize, cTau )

%
% Splits cSignal into pieces
%

if nargin<3,
    cTau = 1;
end;

lN = floor(length(cSignal)/cTau);

vPatches = zeros(lN,cWindowSize);

for k = 1:lN,
    for i = 0:cWindowSize-1,
        idx(i+1) = mod(mod(k-floor(cWindowSize/2),length(cSignal))+i,length(cSignal));
        if idx(i+1)>=0,
            idx(i+1) = idx(i+1)+1;
        end;
    end;
    vPatches(k,:) = cSignal(idx);
end;

return;