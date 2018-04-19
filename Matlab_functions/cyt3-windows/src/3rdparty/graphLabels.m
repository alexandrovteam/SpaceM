function graphLabels( titleLabel, xlabelLabel, ylabelLabel, fontSize, fontWeight )
% graphLabels( titleLabel, xlabelLabel, ylabelLabel, fontSize = 14, fontWeight = 'bold' )
% sets up graph labels according to my favorite style (size 14, bold).

if( nargin < 5 )
	fontWeight = 'bold';
end

if( nargin < 4 )
	fontSize = 14;
end

if( nargin >= 1 )
	title( titleLabel, 'FontSize', fontSize, 'FontWeight', fontWeight )
end

if( nargin >= 2 )
	xlabel( xlabelLabel, 'FontSize', fontSize, 'FontWeight', fontWeight )
end

if( nargin >= 3 )
	ylabel( ylabelLabel, 'FontSize', fontSize, 'FontWeight', fontWeight )
end


