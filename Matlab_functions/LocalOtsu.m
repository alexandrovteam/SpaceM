function oneThresholdedPixel = LocalOtsu(grayImagePatch)
	oneThresholdedPixel = false;
	try
		[rows, columns] = size(grayImagePatch);
		middleRow = ceil(rows/2);
		middleColumn = ceil(columns/2);
		level = graythresh(grayImagePatch);
		% Threshold the center pixel only.
		oneThresholdedPixel = ~im2bw(grayImagePatch(middleRow, middleColumn), level);
	catch ME
		errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
			ME.stack(1).name, ME.stack(1).line, ME.message);
		fprintf(1, '%s\n', errorMessage);
		uiwait(warndlg(errorMessage));
	end
return;