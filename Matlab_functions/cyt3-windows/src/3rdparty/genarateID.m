%Unique id generator 
function uniqueID = genarateID
    rng('shuffle')
    symbolsLetters = 'A':'Z';
    symbolsNumbers = '0':'9';
    nums = randi(numel(symbolsNumbers),[1 5]);
    lets = randi(numel(symbolsLetters),[1 4]);
    nums = symbolsNumbers (nums);
    str = symbolsLetters (lets);
    uniqueID= [str,nums];
end