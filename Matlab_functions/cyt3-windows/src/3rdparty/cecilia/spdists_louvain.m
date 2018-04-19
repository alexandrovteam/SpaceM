function [c, Q] = spdists_louvain( spdists )
% [c, Q] = spdists_louvain( spdists )
%
% Run the C++ version of the Louvain method for community detection. c is a cell array of length d, where d is the number of levels in
% the tree. c{i} is the assignment in level d for each node. Q is the modularity at the end of each iteration.
%
% spdists is expected to be undirected. The script creates several temporary files which are not removed afterwards.
%
% NOTE: the weight of the edges is irrelevant (they are always exported as weight 1!).

[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
curr_path = [curr_path filesep]

% generate a random file name to use using current time
filename = [curr_path 'spdists_' num2str( now ) ];


% export spdists to a text file
fprintf( 1, 'MATLAB: exporting %s.txt\n', filename );

out = fopen( [ filename '.txt' ], 'w' );

f = find( spdists > 0 );
[i, j] = ind2sub( size( spdists ), f );
for idx = 1:length( i )
	fprintf( out, '%d %d\n', i( idx )-1, j( idx )-1 );
end

fclose( out );


% convert txt file to bin file
fprintf( 1, 'MATLAB: calling convert:\n' );
command = [curr_path 'convert -i ' filename '.txt -o ' filename '.bin' ];
fprintf( 1, '%s\n', command );
system( command );
fprintf( 1, '\n' );


% run louvain algorithm
fprintf( 1, 'MATLAB: running louvain algorithm:\n' );
command = [curr_path 'community ' filename '.bin -l -1 -v > ' filename '.graph.tree' ];
fprintf( 1, '%s\n', command );
[s, r] = system( command );
fprintf( 1, '\n' );
% find each iteration's modularity
fprintf( 1, 'MATLAB: modularity scores:\n' );
Q = find_modularity( r )

% find nu. of levels
command = [curr_path 'hierarchy ' filename '.graph.tree' ];
fprintf( 1, '%s\n', command );
[s, r] = system( command );
fprintf( 1, '\n' );

r = strtok( r, 10 );
r = regexprep( r, 'Number of levels: ', '' );
nu_levels = str2num( r )-1;


fprintf( 1, 'MATLAB: max level is %d\n', nu_levels );


% import each of the levels g.t. 0
for level = 1:nu_levels
	fprintf( 1, 'MATLAB: importing level %d\n', level );

	command = [curr_path 'hierarchy ' filename '.graph.tree -l ' num2str( level ) ' > ' filename '.tmp' ];
	fprintf( 1, '%s\n', command );
	system( command );

	hierarchy_output = load( [ filename '.tmp' ] );

	c{level} = hierarchy_output( :, 2 ) + 1;
end


	function Q = find_modularity( r )
	% Q = find_modularity( r )
	%
	% convert the text output into modularity score of each iteration
	signature = '  modularity increased from %f to %f';
	idx = 0;

	while( ~isempty( r ) )
		% read a line and match it to the signature
		[token, r] = strtok( r, char( 10 ) );
		a = sscanf( token, signature );

		if( ~isempty( a ) )
			% signature matched copy new modularity
			idx = idx + 1;
			Q( idx ) = a( 2 );
		end
	end


