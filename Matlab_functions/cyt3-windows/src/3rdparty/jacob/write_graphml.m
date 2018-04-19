function write_graphml(spdists, sampleID, samplePercent, means, names, cluster_names )
% spdists_export_graphml( filename, spdists )
%
% export spdists as a graphml file

% get filename to save to
[filename,pathname,~] = uiputfile('*.graphml','Save Graph');

if isequal(filename,0) || isequal(pathname,0)
    return;
end
filename = [pathname filename];

spdists = triu( spdists ); % prevent double edges by examining triu
n = length( spdists );

out = fopen( filename, 'w' );

% header
fprintf( out, '<?xml version="1.0" encoding="UTF-8"?>\n<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">\n' );
fprintf( out, '<key id="weight" for="edge" attr.name="weight" attr.type="double"/>\n' );
%node attributes
fprintf( out, '<key id="sampleID" for="node" attr.name="sampleID" attr.type="double"/>\n' );
fprintf( out, '<key id="samplePercent" for="node" attr.name="samplePercent" attr.type="double"/>\n' );
for j = 1:size(means,2)
    fprintf( out, '<key id="%s" for="node" attr.name="%s" attr.type="double"/>\n', names{j}, names{j} );
end
fprintf( out, '<graph id="G" edgedefault="undirected">\n' );

% export nodes
for i = 1:n
	fprintf( out, '\t<node id="n%d">\n', i );
	fprintf( out, '\t<node label="%s">\n', cluster_names{i} );
    fprintf(out, '\t<data key="sampleID">%i</data>\n', sampleID(i) );
    fprintf(out, '\t<data key="samplePercent">%i</data>\n', samplePercent(i) );
    for j = 1:size(means,2)
        fprintf(out, '\t<data key="%s">%3.2f</data>\n', names{j}, means(i,j) );
    end
	fprintf( out, '\t</node>\n' );
end

fprintf( out, '\n' );

% export edges
[ i, j, s ] = find( spdists > 0 );

for idx = 1:length( i )
	fprintf( out, '\t<edge id="e%d" source="n%d" target="n%d">\n', idx, i( idx ), j( idx ) );
	fprintf( out, '\t\t<data key="weight">%3.4f</data>\n', s( idx ) );
	fprintf( out, '\t</edge>\n' );

	if( mod( idx, 10000 ) == 0 )
		fprintf( 1, '%3.2f%%\n', idx / length( i ) * 100 );
	elseif( mod( idx, 500 ) == 0 )
		fprintf( 1, '.' );
	end
end

fprintf( out, '</graph>\n</graphml>\n' );

fclose( out );