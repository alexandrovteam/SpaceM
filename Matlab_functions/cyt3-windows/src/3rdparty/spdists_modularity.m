function Q = spdists_modularity( spdists, c, verbose )
% Q = spdists_modularity( spdists, c, verbose )
%
% calculate modularity of spdists given partition c. verbose is 'off' by default.
%
% calculation is done iteratively (not through matrix calculation), therefore it is quite slow.

if( ~exist( 'verbose', 'var' ) )
	verbose = 0;
end

if( strcmpi( verbose, 'on' ) )
	verbose = 1;
end

% calculate degree of each node
deg = sum( spdists );
m2 = sum( deg ); % 2 * m

% find indices of each community
for comm_idx = 1:max( c )
	comm_indices{comm_idx} = find( c == comm_idx );
end

% calculate A_ij - ( k_i k_j ) / 2m, for each pair of nodes in the same community
s = 0;
for comm_idx = 1:max( c )
	if( verbose )
		fprintf( 1, 'calculating sum for community n%d (%d nodes)\n', comm_idx, length( comm_indices{comm_idx} ) );
	end

	t = tic;

	s = s + sum( sum( spdists( comm_indices{comm_idx}, comm_indices{comm_idx} ) ) ) - v_prime_times_v( deg( comm_indices{comm_idx} ), verbose ) / m2;

	%{
		for i_idx = 1:length( comm_indices{comm_idx} )
			i = comm_indices{comm_idx}( i_idx );
			for j_idx = i_idx+1:length( comm_indices{comm_idx} )
				j = comm_indices{comm_idx}( j_idx );

				s = s + spdists( i, j ) - deg( i ) * deg( j ) / m2;
			end
		end
	%}

	toc( t );
	fprintf( 1, '\n' );
end

Q = s / m2;


	function s = v_prime_times_v( v, verbose )
	% s = v_prime_times_v( v, verbose )
	%
	% calculate sum( v' * v ).
	s = 0;

	for i = 1:length( v )
		s = s + sum( v( i ) * v );

		% print user status report
		if( verbose )
			if( rem( i, 1000 ) == 0 )
				fprintf( 1, '%3.2f%%\n', i / length( v ) * 100 );
			elseif( rem( i, 25 ) == 0 )
				fprintf( 1, '.' );
			end
		end
	end


