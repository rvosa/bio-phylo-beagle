package Bio::Phylo::Beagle;
use strict;
use warnings;
use beagle;
use Bio::Phylo::BeagleOperation;
use Bio::Phylo::BeagleOperationArray;
use Data::Dumper;
use base 'Bio::Phylo';
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Phylo::Util::CONSTANT qw':objecttypes /looks_like/';
use Bio::Phylo::Util::Logger;

my $logger = Bio::Phylo::Util::Logger->new;
our $BEAGLE_OP_NONE = -1;
my ( %model, %matrix, %tree, %instance );

sub new {
    $logger->info("@_");
	my $self = shift->SUPER::new(@_);
	#unless( $self->get_tree && $self->get_matrix && $self->get_model ) {
	#	throw 'BadArgs' => 'Need tree, matrix and model';
	#}
	return $self;
}

# static
my $create_intarray = sub {
    $logger->debug("create_intarray: @_");
	my @ints = @_;
	my $result = beagle::new_intArray(scalar(@ints));
	for my $i ( 0 .. $#ints ) {
		beagle::intArray_setitem($result,$i,$ints[$i]);
	}
	return $result;
};

# static
my $create_doublearray = sub {
    $logger->debug("create_doublearray: @_");
	my @doubles = @_;
	my $result = beagle::new_doubleArray(scalar(@doubles));
	for my $i ( 0 .. $#doubles ) {
		beagle::doubleArray_setitem($result,$i,$doubles[$i]);
	}
	return $result;
};

# static
my $create_states = sub {
    $logger->debug("create_states: @_");
	my ( $seq, $table ) = @_;
	my @char = $seq->get_char;
	my $states = beagle::new_intArray(scalar(@char));
	for my $i ( 0 .. $#char ) {
		beagle::intArray_setitem($states,$i,$table->{$char[$i]});
	}
	return $states;
};

# static
my $create_pattern_weights = sub {
    $logger->debug("create_pattern_weights: @_");
	my @wts = @_;
	my $weights = beagle::new_doubleArray(scalar(@wts));
	for my $i ( 0 .. $#wts ) {
		beagle::doubleArray_setitem($weights,$i,$wts[$i]);
	}
	return $weights;
};

# method
sub create_table {
    $logger->info("@_");
	my $self   = shift;
	my $matrix = $self->get_matrix;
    my %seen;
    for my $row ( @{ $matrix->get_entities } ) {
        my @char = $row->get_char;
        $seen{$_}++ for @char;
    }
    my @states = sort { $a cmp $b } grep { /[A-Z]/ } keys %seen;
    my %table;
    my $counter = 0;
    $table{$_} = $counter++ for @states;
    for my $key ( keys %table ) {
        $table{lc $key} = $table{$key};
    }
    $table{'-'} = scalar(@states);
    return \%table;
}

=item set_pattern_weights

 Type    : Mutator
 Title   : set_pattern_weights
 Usage   : $beagle->set_pattern_weights( 1,1,1,1,1 )
 Function: Set a category weights buffer
 Returns : error code
 Args    : Array containing patternCount weights (input)

=cut

sub set_pattern_weights {
    $logger->info("@_");
	my $self     = shift;
	my $nchar    = $self->get_matrix->get_nchar;
	my $instance = $self->get_instance;
    my @args = scalar(@_) ? @_ : ( 1 );
	my @wts;
    for my $i ( 1 .. $nchar ) {
        push @wts, $args[$i] || 1;
    }
	my $patternWeights = $create_pattern_weights->(@wts);
    
    #/**
    #* @brief Set pattern weights
    #*
    #* This function sets the vector of pattern weights for an instance.
    #*
    #* @param instance              Instance number (input)
    #* @param inPatternWeights      Array containing patternCount weights (input)
    #*
    #* @return error code
    #*/
	return beagle::beagleSetPatternWeights($instance, $patternWeights);
}

=item set_state_frequencies

This function copies a state frequency array into an instance buffer.

 Type    : Mutator
 Title   : set_state_frequencies
 Usage   : $beagle->set_state_frequencies
 Function: Set a state frequency buffer
 Returns : error code
 Args    : Optional: Index of state frequencies buffer (input)

=cut

sub set_state_frequencies {
    $logger->info("@_");
	my $self  = shift;
    my $index = shift || 0;
	my $model = $self->get_model;
	my $inst  = $self->get_instance;
	my $table = $self->create_table;
	my %keys  = map { uc $_ => 1 } keys %{ $table };
	delete $keys{'-'};
	delete $keys{'?'};
	my @states = sort { $a cmp $b } keys %keys;
	my @base_freqs;
	push @base_freqs, $model->get_pi($_) for @states;
	my $freqs = $create_pattern_weights->(@base_freqs);
    
    #/**
    # * @brief Set a state frequency buffer
    # *
    # * This function copies a state frequency array into an instance buffer.
    # *
    # * @param instance              Instance number (input)
    # * @param stateFrequenciesIndex Index of state frequencies buffer (input)
    # * @param inStateFrequencies    State frequencies array (stateCount) (input)
    # *
    # * @return error code
    # */    
	return beagle::beagleSetStateFrequencies($inst,$index,$freqs);
}

=item create_instance

This function creates a single instance of the BEAGLE library and can be
called multiple times to create multiple data partition instances each
returning a unique identifier.

 Type    : Mutator
 Title   : create_instance
 Usage   : $beagle->create_instance( %args )
 Function: Create a single instance
 Returns : the unique instance identifier (<0 if failed, see @ref BEAGLE_RETURN_CODES "BeagleReturnCodes")
 Args    : -tip_count              => Number of tip data elements (input)
           -partials_buffer_count  => Number of partials buffers to create (input)
           -compact_buffer_count   => Number of compact state representation buffers to create (input)
           -state_count            => Number of states in the continuous-time Markov chain (input)
           -pattern_count          => Number of site patterns to be handled by the instance (input)
           -eigen_buffer_count     => Number of rate matrix eigen-decomposition, category weight, and state frequency buffers to allocate (input)
           -matrix_buffer_count    => Number of transition probability matrix buffers (input)
           -category_count         => Number of rate categories (input)
           -scale_buffer_count		=> Number of scale buffers to create, ignored for auto scale or always scale (input)
           -resource_list          => List of potential resources on which this instance is allowed (input, NULL implies no restriction)
           -resource_count         => Length of resourceList list (input)
           -preference_flags       => Bit-flags indicating preferred implementation characteristics, see BeagleFlags (input)
           -requirement_flags      => Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
           -return_info            => Pointer to return implementation and resource details

=cut

sub create_instance {
    $logger->info("@_");
    my $self = shift;
    if ( my %args = looks_like_hash @_ ) {
        my $matrix = $self->get_matrix;
        my $info = beagle::BeagleInstanceDetails->new();
        my $instance = beagle::beagleCreateInstance(
            $args{'-tip_count'}             || $matrix->get_ntax, # tipCount
            $args{'-partials_buffer_count'} || 0, # partialsBufferCount
            $args{'-compact_buffer_count'}  || 0, # compactBufferCount
            $args{'-state_count'}           || 0, # stateCount
            $args{'-pattern_count'}         || $matrix->get_nchar, # patternCount
            $args{'-eigen_buffer_count'}    || 0, # eigenBufferCount
            $args{'-matrix_buffer_count'}   || 0, # matrixBufferCount
            $args{'-category_count'}        || 0, # categoryCount
            $args{'-scale_buffer_count'}    || 0, # scaleBufferCount
            $args{'-resource_list'}         || undef, # resourceList
            $args{'-resource_count'}        || 0, # resourceCount
            $args{'-preference_flags'}      || 0, # preferenceFlags,
            $args{'-requirements_flags'}    || 0, # requirementsFlags,
            $args{'-return_info'}           || $info # returnInfo
        );
        $instance{$self->get_id} = $instance;
        if ( $instance < 0 ) {
            throw 'ExtensionError' => "Can't instantiate BEAGLE";
        }
        return $instance;
    }
}

sub get_instance { $instance{ shift->get_id } }

sub get_tree { $tree{ shift->get_id } }

sub set_tree {
    $logger->info("@_");
	my ( $self, $tree ) = @_;
	if ( looks_like_object $tree, _TREE_ ) {
		$tree{ $self->get_id } = $tree;
		return $self;
	}
}

sub get_matrix { $matrix{ shift->get_id } }

sub set_matrix {
    $logger->info("@_");
	my ( $self, $matrix ) = @_;
	if ( looks_like_object $matrix, _MATRIX_ ) {
		$matrix{ $self->get_id } = $matrix;
		my $table = $self->create_table;
		my $instance = $self->get_instance;
		my @seqs = @{ $matrix->get_entities };
		for my $i ( 0 .. $#seqs ) {
			my $states = $create_states->($seqs[$i],$table);
			beagle::beagleSetTipStates($instance,$i,$states);
		}
		$self->set_pattern_weights;
		return $self;
	}
}

sub get_model { $model{ shift->get_id } }

sub set_model {
    $logger->info("@_");
	my ( $self, $model ) = @_;
	if ( looks_like_object $model, _MODEL_ ) {
		$model{ $self->get_id } = $model;
		$self->set_state_frequencies;
		return $self;
	}
}

=item set_category_weights

This function copies a category weights array into an instance buffer.

 Type    : Mutator
 Title   : set_category_weights
 Usage   : $beagle->set_category_weights( -weights => [1.0] )
 Function: Set a category weights buffer
 Returns : error code
 Args    : -weights => [ Category weights array (categoryCount) (input) ]
           -index => Optional: index of category weights buffer

=cut

sub set_category_weights {
    $logger->info("@_");
    my $self = shift;
    if ( my %args = looks_like_hash @_ ) {
        
        # /**
        #  * @brief Set a category weights buffer
        #  *
        #  * This function copies a category weights array into an instance buffer.
        #  *
        #  * @param instance              Instance number (input)
        #  * @param categoryWeightsIndex  Index of category weights buffer (input)
        #  * @param inCategoryWeights     Category weights array (categoryCount) (input)
        #  *
        #  * @return error code
        #  */
        my $weights  = $create_pattern_weights->(@{ $args{'-weights'} });
        my $instance = $self->get_instance;
        my $index = $args{'-index'} || 0;
        return beagle::beagleSetCategoryWeights($instance, $index, $weights);
    }
}

=item set_category_rates

This function sets the vector of category rates for an instance.

 Type    : Mutator
 Title   : set_category_rates
 Usage   : $beagle->set_category_rates( 1.0 )
 Function: Set category rates
 Returns : error code
 Args    : Array containing categoryCount rate scalers (input)

=cut

sub set_category_rates {
    $logger->info("@_");
    my ( $self, @args ) = @_;
    my @rates = scalar(@args) || (1.0);
    # /**
    #  * @brief Set category rates
    #  *
    #  * This function sets the vector of category rates for an instance.
    #  *
    #  * @param instance              Instance number (input)
    #  * @param inCategoryRates       Array containing categoryCount rate scalers (input)
    #  *
    #  * @return error code
    #  */
    my $rates = $create_pattern_weights->(@rates);
    my $instance = $self->get_instance;
    return beagle::beagleSetCategoryRates($instance, $rates);
}

=item set_eigen_decomposition

This function copies an eigen-decomposition into an instance buffer.

 Type    : Mutator
 Title   : set_eigen_decomposition
 Usage   : $beagle->set_eigen_decomposition( %args )
 Function: Set an eigen-decomposition buffer
 Returns : error code
 Args    : -index           => Optional: index of eigen-decomposition buffer           
           -vectors         => Flattened matrix (stateCount x stateCount) of eigen-vectors (input)
           -inverse_vectors => Flattened matrix (stateCount x stateCount) of inverse-eigen- vectors (input)
           -values          => Vector of eigenvalues

=cut

sub set_eigen_decomposition {
    $logger->info("@_");
    my $self = shift;
    if ( my %args = looks_like_hash @_ ) {
        my $evec  = $create_pattern_weights->( @{ $args{'-vectors'} } );
        my $ivec  = $create_pattern_weights->( @{ $args{'-inverse_vectors'} } );
        my $eval  = $create_pattern_weights->( @{ $args{'-values'} } );
        my $index = $args{'-index'} || 0;
        my $instance = $self->get_instance;
        
        # /**
        #  * @brief Set an eigen-decomposition buffer
        #  *
        #  * This function copies an eigen-decomposition into an instance buffer.
        #  *
        #  * @param instance              Instance number (input)
        #  * @param eigenIndex            Index of eigen-decomposition buffer (input)
        #  * @param inEigenVectors        Flattened matrix (stateCount x stateCount) of eigen-vectors (input)
        #  * @param inInverseEigenVectors Flattened matrix (stateCount x stateCount) of inverse-eigen- vectors
        #  *                               (input)
        #  * @param inEigenValues         Vector of eigenvalues
        #  *
        #  * @return error code
        #  */
        return beagle::beagleSetEigenDecomposition($instance,$index,$evec,$ivec,$eval);
    }
}

=item update_transition_matrices

This function calculates a list of transition probabilities matrices and their
first and second derivatives (if requested).

 Type    : Mutator
 Title   : update_transition_matrices
 Usage   : $beagle->update_transition_matrices( %args )
 Function: Calculate a list of transition probability matrices
 Returns : error code
 Args    : -index  => Optional: Index of eigen-decomposition buffer
           -deriv1 => Optional: List of indices of first derivative matrices to update
           -deriv2 => Optional: List of indices of second derivative matrices to update

=cut

sub update_transition_matrices {
    $logger->info("@_");
    my $self = shift;
    my %args = looks_like_hash @_;
        
    # create node and edge arrays
    my $tree = $self->get_tree;        
    my ( @nodeIndices, @edgeLengths );
    my $nodeIndex = 0;
    $tree->visit_depth_first(
        '-post' => sub {
            my $node = shift;
            if ( not $node->is_root ) {
                push @nodeIndices, $nodeIndex;
                push @edgeLengths, $node->get_branch_length;
                $nodeIndex++;
            }
        }
    );
    my $nodeIndices = $create_intarray->(@nodeIndices);
    my $edgeLengths = $create_doublearray->(@edgeLengths);
    
    my $instance = $self->get_instance;
    my $index = $args{'-index'} || 0;
    my $firstderiv  = $args{'-deriv1'}; # XXX maybe do with $create_intarray
    my $secondderiv = $args{'-deriv2'};
    
    # /**
    #  * @brief Calculate a list of transition probability matrices
    #  *
    #  * This function calculates a list of transition probabilities matrices and their first and
    #  * second derivatives (if requested).
    #  *
    #  * @param instance                  Instance number (input)
    #  * @param eigenIndex                Index of eigen-decomposition buffer (input)
    #  * @param probabilityIndices        List of indices of transition probability matrices to update
    #  *                                   (input)
    #  * @param firstDerivativeIndices    List of indices of first derivative matrices to update
    #  *                                   (input, NULL implies no calculation)
    #  * @param secondDerivativeIndices    List of indices of second derivative matrices to update
    #  *                                   (input, NULL implies no calculation)
    #  * @param edgeLengths               List of edge lengths with which to perform calculations (input)
    #  * @param count                     Length of lists
    #  *
    #  * @return error code
    #  */        
    return beagle::beagleUpdateTransitionMatrices(
        $instance,    # instance
        0,            # eigenIndex
        $nodeIndices, # probabilityIndices
        $firstderiv,  # firstDerivativeIndices
        $secondderiv, # secondDerivativeIndices
        $edgeLengths, # edgeLengths
        scalar(@nodeIndices) # count
    );
}

#my $operations = Bio::Phylo::BeagleOperationArray->new(2);
#
#my $op0 = Bio::Phylo::BeagleOperation->new(
#    '-destination_partials'     => 3,
#    '-destination_scale_write'  => $BEAGLE_OP_NONE,
#    '-destination_scale_read'   => $BEAGLE_OP_NONE,
#    '-child1_partials'          => 0,
#    '-child1_transition_matrix' => 0,
#    '-child2_partials'          => 1,
#    '-child2_transition_matrix' => 1
#);
#
#my $op1 = Bio::Phylo::BeagleOperation->new(
#    '-destination_partials'     => 4,
#    '-destination_scale_write'  => $BEAGLE_OP_NONE,
#    '-destination_scale_read'   => $BEAGLE_OP_NONE,
#    '-child1_partials'          => 2,
#    '-child1_transition_matrix' => 2,
#    '-child2_partials'          => 3,
#    '-child2_transition_matrix' => 3
#);
#$operations->set_item( -index => 0, -op => $op0 );
#$operations->set_item( -index => 1, -op => $op1 );

sub update_partials {
    $logger->info("@_");
    my $self = shift;
    if ( my %args = looks_like_hash @_ ) {
        
        my $operations = $args{'-operations'} || throw 'BadArgs' => 'Need -operations argument';
        my $count      = $args{'-count'}      || throw 'BadArgs' => 'Need -count argument';
        my $index      = $args{'-index'}      || throw 'BadArgs' => 'Need -index argument';
        
        # /**
        #  * @brief Calculate or queue for calculation partials using a list of operations
        #  *
        #  * This function either calculates or queues for calculation a list partials. Implementations
        #  * supporting ASYNCH may queue these calculations while other implementations perform these
        #  * operations immediately and in order.
        #  *
        #  * @param instance             Instance number (input)
        #  * @param operations           BeagleOperation list specifying operations (input)
        #  * @param operationCount       Number of operations (input)
        #  * @param cumulativeScaleIndex Index number of scaleBuffer to store accumulated factors (input)
        #  *
        #  * @return error code
        #  */
        beagle::beagleUpdatePartials(
            $self->get_instance,
            $operations->get_array,
            $count,
            $index
        );
    }
}

sub calculate_root_log_likelihoods {
    $logger->info("@_");
    my $self = shift;
    my %args = looks_like_hash @_;
    my $outSumLogLikelihood = beagle::new_doublep();
    my $categoryWeightsIndices = $create_intarray->($args{'-category_weights_indices'} || 0);
    my $stateFrequencyIndices  = $create_intarray->($args{'-state_frequencies_indices'} ||0);
    my $cumulativeScaleIndices = $create_intarray->($BEAGLE_OP_NONE);
    
    #/**
    # * @brief Calculate site log likelihoods at a root node
    # *
    # * This function integrates a list of partials at a node with respect to a set of partials-weights
    # * and state frequencies to return the log likelihood sum
    # *
    # * @param instance                 Instance number (input)
    # * @param bufferIndices            List of partialsBuffer indices to integrate (input)
    # * @param categoryWeightsIndices   List of weights to apply to each partialsBuffer (input). There
    # *                                  should be one categoryCount sized set for each of
    # *                                  parentBufferIndices
    # * @param stateFrequenciesIndices  List of state frequencies for each partialsBuffer (input). There
    # *                                  should be one set for each of parentBufferIndices
    # * @param cumulativeScaleIndices   List of scaleBuffers containing accumulated factors to apply to
    # *                                  each partialsBuffer (input). There should be one index for each
    # *                                  of parentBufferIndices
    # * @param count                    Number of partialsBuffer to integrate (input)
    # * @param outSumLogLikelihood      Pointer to destination for resulting log likelihood (output)
    # *
    # * @return error code
    # */
    my $bufferIndices = $create_intarray->( scalar( @{ $self->get_tree->get_entities } ) - 1 );
    beagle::beagleCalculateRootLogLikelihoods(
        $self->get_instance,
        $bufferIndices,
        $categoryWeightsIndices,
        $stateFrequencyIndices,
        $cumulativeScaleIndices,
        $args{'-count'} || 1,
        $outSumLogLikelihood
    );
    # 
    return beagle::doublep_value($outSumLogLikelihood);    
}

1;