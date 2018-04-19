function spdists = spdists_knngraph( x, k, distance, chunk_size, verbose )
% spdists = spdists_knngraph( x, k, distance, chunk_size, verbose )
%
% finds the nearest neighbor graph for data matrix x. the graph is represented as a sparse matrix. for each point in x,
% the k nearest neighbors are found. distance is the distance metric used by knnsearch (such as euclidean, cosine,
% etc.).
%
% if chunk_size is specified, x will be split to chunk_size-by-D submatrices, and the calculation will be run
% separately on each submatrix (in order to conserve memory).

	if( ~exist( 'distance', 'var' ) )
		distance = 'euclidean';
	end

	if( ~exist( 'chunk_size', 'var' ) )
		chunk_size = size( x, 1 );
	end

	if( ~exist( 'verbose', 'var' ) )
		verbose = 0;
	end

	n = size( x, 1 );

	spdists = sparse( n, n );

	chunk_num = 0;
	total_chunks = ceil( n / chunk_size );

	% iterate over submatrices
	for from = 1:chunk_size:n
		to = min( from + chunk_size - 1, n );
		rx = from:to;
		x_from_to = x( rx, : );
		chunk_num = chunk_num + 1;

		[ idx, d ] = knnsearch( x, x_from_to, 'k', k + 1, 'distance', distance );

		idx( :, 1 ) = []; d( :, 1 ) = []; % remove self neighbor

		% update spdists
		js = repmat( rx', 1, k );
		indices = sub2ind( size( spdists ), idx(:), js(:) );
		spdists( indices ) = d(:);

		if( verbose )
			% report progress if required by user
			if( mod( chunk_num, 10 ) == 0 )
				fprintf( 1, '.' );
			elseif( mod( chunk_num, 100 ) == 0 )
				fprintf( 1, '%3.2f%%', chunk_num / total_chunks * 100 );
			end
		end
	end


