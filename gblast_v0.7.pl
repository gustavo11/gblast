#!/bin/env perl

use strict;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Getopt::Std;
#use Term::ProgressBar;

my $MAX_E_VALUE = 1000;
my $DEFAULT_WIDTH = 1000;
my $DEFAULT_LENGTH = 1500;
my $DEFAULT_COLOR = "black";
my $LEFT_BORDER = 0;
STDOUT->autoflush(1);


my $usage = "\n\nv0.7\ngblast -q <FASTA of query sequences>  -r <m8 blast file> [ -s <FASTA of subject sequences> ] [-o <output png file>] [-w <width>] [-e <e-value threshold> ] [-c <coord. tolerance>] [ -g <group file>]\n\n" .
"-o output file. default output name = <input file name>_[<group info>|<e_value>]_<coord. tolerance>.svg\n" .
"-w width of final figure. default=$DEFAULT_WIDTH.\n" .
"-s FASTA of subject sequences. If not provided the figure will not show the length of subject sequence associated with each HSP.\n" .
"-e quick way to set an overall cutoff for e-value. The e-value thresholds in the \"group file\" have precedence over this parameter.\n" .
"-c <coord. tolerance>\t. Indicates that HSP with higher scores should have priority during rendering." .
"HSPs of lower score will not be rendered if they are <coord. tolerance> nt apart from HSPs with higher score. distance = max( start - start', end - end' )\n" .
"Obs.: \'-c 0\' is different than omiting this parameter. \'-c 0\' means that HSPs of lower score will not be rendered if they overlap (distance = 0) with HSP of higher scores\n" .
"If the '-c' parameter is omited, HSP\'s covering the exact same region will be shown.\n\n" .  
"-g <group file>\tdefines group of subject sequences based on regex. E-value threshold and glyph color can be assigned to each group. \n" .
"The definition file is a tab delimited file using the following format:\n\n" .
"<group name>\t[\\w\\W]+\t\t\t<color ex.:red,green,blue>\t<e-value threshold ex.:0.1, 0.000001>\n" .
"<group name>\tgi\\|74584587\\|809998\t<color ex.:red,green,blue>\t<e-value threshold ex.:0.1, 0.000001>\n\n" .
"Obs.: Please escape metacharacters in the regex, such as the \'|\' in the previous example\n\n" .
"Output description:\n" .
"-------------------\n" .
"Description of items in the red caption associated to each HSP:\n".
"qc = coordinates in the query sequence.\n" .
"sc = coordinates in the subject sequence.\n" .
"sl = subject length.\n" .
"s  = score.\n" .
"l  = length of alignment in the query sequence.\n" .
"e  = e-value.\n" .
"%i = percentage of identity of the aligned fragment.\n\n\n";

##############################################################################
# Dealing with command line parameters
##############################################################################
my %opt;

#getopts('oif:', \%opts);  o e i sao flags e f eh argumento

getopts(':q:w:r:o:e:c:g:s:',\%opt);

die $usage if( ! defined( $opt{r} ) );
die $usage if( ! defined( $opt{q} ) );
die $usage if( ! defined( $opt{s} ) );

my $width             = $opt{w};
my $queryFasta        = $opt{q};
my $subjectFasta      = $opt{s};
my $inputFile         = $opt{r};
my $outputFile        = $opt{o};
my $eThreForAll       = $opt{e};
my $coordTolerance    = $opt{c};
my $grpFile           = $opt{g};
my $usePages          = $opt{p};

my @group;

# Caso o usuario nao esteje utilizando o argumento
# para definica de valor de corte de e-value
$eThreForAll = $MAX_E_VALUE if( !defined $eThreForAll );

# Caso o usuario esteje utilizando um arquivo
# de definicao de grupos
if(  defined $grpFile ){
  print "Using \'$grpFile\' as group definition...\n";
  open GRP, $grpFile;

  my $ind = 0;
  while( <GRP> ){
    my $line = $_;
    if( my ($name,$regexp,$color,$e_threshold) = split /[\t\s]+/ ){
      
      $group[$ind]{name} 	= $name;
      $group[$ind]{regexp}      = $regexp;
      $group[$ind]{color}       = $color;
      $group[$ind]{e_threshold} = $e_threshold;
      $ind++;
    }else{
      NOT_RECON_GRP:
      $line =~ s/[\n\r]//g;
      die "Not recognized group file \'$grpFile\'. Error on line " . ($ind + 1) . " : \'$line\'\n";
    }
  } 

# Caso o usuario nao esteje utilizando um arquivo
# de grupos
}else{
  $group[0]{color} = "blue";
  $group[0]{regexp} = "[\\w\\W]+";
  $group[0]{e_threshold} = $eThreForAll;
}

$width = $DEFAULT_WIDTH if( !defined $width );

############################################################################

# Count the number of lines in the file
my ($empty,$lines,$words,$letters) = 
  split /\s+/, qx{wc $inputFile};

  
# Read FASTA files

# Query
my %query_length;
my $inSeqIO     = Bio::SeqIO->new(-file => $queryFasta, '-format' => 'Fasta');
while ( my $inSeq = $inSeqIO->next_seq() ){
  if( !defined( $inSeq->seq() ) || length( $inSeq->seq() ) == 0 ){
  	  print STDERR  $inSeq->id() . " have an empty string as sequence. Sequence discarded\n";
  } else {
  	  $query_length{ $inSeq->id() } = $inSeq->length();
  }
}
$inSeqIO->close();

# Subject
my %subject_length;
if( defined $subjectFasta ){
my $inSeqIO     = Bio::SeqIO->new(-file => $subjectFasta, '-format' => 'Fasta');
while ( my $inSeq = $inSeqIO->next_seq() ){
  if( !defined( $inSeq->seq() ) || length( $inSeq->seq() ) == 0 ){
  	  print STDERR  $inSeq->id() . " have an empty string as sequence. Sequence discarded\n";
  } else {
  	  $subject_length{ $inSeq->id() } = $inSeq->length();
  }
}
$inSeqIO->close();
}


  
# Retrieve the length of the largest query sequence
my $maxLength = 0;
#my $progMaxLength = Term::ProgressBar->new($lines);
my $contIteMaxLength = 0;
open INPUT_MAX_LENGTH, $inputFile;
while (<INPUT_MAX_LENGTH>) {		 
  my $line = $_;

  chomp $line;
  next if $line =~ /^\#/;  		 # ignore comments
  my( $query, $subject, $perc_ident, $reg_sim_length, $num_mismatches, $num_gaps, $q_start, $q_end, $s_start, $s_end, $e_value, $score ) = split "\t", $line;
  
  die "Unable to find the query sequence $query in the FASTA file!\n" if ( ! defined( $query_length{ $query } ) );

  $maxLength = $query_length{ $query } if( $query_length{ $query } > $maxLength );
  $contIteMaxLength++;
#  $progMaxLength->update( $contIteMaxLength );
}
#$progMaxLength->update($lines);
close INPUT_MAX_LENGTH;

$maxLength += $LEFT_BORDER;

my $panel = Bio::Graphics::Panel->new(-length     => $maxLength,
				      -image_class=>'GD::SVG', 
				      -width      => $width,
				      -pad_left   => 40,
				      -pad_right  => 200,
				      -pad_top    => 40,
				      -pad_bottom => 40,
				     );

my $old_name_query = "";
my $track;

print "Sorting by decreasing score value...\n";
  
my @blast_result;  
open INPUT, $inputFile;
while( <INPUT> ){
	my $line = $_;
	chomp $line;
	push @blast_result, $line; 	
}
close(INPUT);

# Sort based on query and then based on the score
my @sorted_result = sort { (split /[\t\s]+/, $a)[0] cmp (split /[\t\s]+/, $b)[0] ||  (split /[\t\s]+/, $b)[11] cmp (split /[\t\s]+/, $a)[11] } @blast_result;
  

# Store the coordinates in the query that each of the defined group
# has already covered

# Arrays of Arrays
# first index: group
# second index: coord.
my @qryCov;

print "Drawing graph...\n";

# Progress Counter
#my $progSubjIte = Term::ProgressBar->new($lines);
my $contIteSubj = 0;

SUBJECT_ITERATION:
foreach my $line ( @sorted_result ){
  $contIteSubj++;
#  $progSubjIte->update( $contIteSubj );

 
  my( $query, $subject, $perc_ident, $reg_sim_length, $num_mismatches, $num_gaps, $q_start, $q_end, $s_start, $s_end, $e_value, $score ) = split /[\t\s]+/, $line;
  
  die "Unable to find query sequence $query in the FASTA file!\n" if ( ! defined( $query_length{ $query } ) );
  if( defined $subjectFasta ){
  	die "Unable to find subject sequence $subject in the FASTA file!\n" if ( ! defined( $subject_length{ $subject } ) );
  }

  my $q_length = $query_length{ $query };
  my $s_length = 'NA';
  $s_length = $subject_length{ $subject } if defined $subjectFasta;
  
  #my( $query, $date, $q_length, $blast_v, $q_file, $subject, $q_start, $q_end, $s_start, $s_end,$perc_ident, $unknown1, $score, $unknown2, $unknown3, $s_description, $unknown4, $plus_minus,$s_length, $e_value, $unknown5 ) = split /\//;

  #Add 10 on numbers using scientific format without it (Ex.: e-11 to 10e-11) 
  $e_value =~ s/e/10e/ if( $e_value =~ /^e/ );

  if ( $old_name_query ne $query ) {
    $old_name_query = $query;

    # Cleaning
    undef @qryCov;

    my $full_length = Bio::SeqFeature::Generic->new(-start=>1,-end=>$q_length);

    $panel->add_track($full_length,
		      -glyph   => 'arrow',
		      -tick    => 2,
		      -fgcolor => 'black',
		      -font2color => 'black',
		      -double  => 1,
		      -label  => $query . " size:" . $q_length
		     );

    $track = $panel->add_track(-glyph => 'graded_segments',
				  -label  => 1,

# sub{
#					my $feature = shift;
#                                    	my $name   = $feature->display_name;
#					return $name . " lenght: " . $subject;
#				},
						
				  -bgcolor =>sub {
				    my $feature = shift;
				    my $name = $feature->display_name;

				    my $groupNum;
				    if( ( $groupNum = getGroupNum( $name ) ) != -1 ){
				      return $group[ $groupNum ]->{color};
				    }else{
				      return $DEFAULT_COLOR ;
				    }
	  
				   },
				  -min_score => 0,
				  -max_score => 0,
				  -font2color     => 'red',
				  -sort_order     => 'high_score',
				  -description => sub {
				    my $feature = shift;
				    #my $sc   = $feature->score;
				    my $e    = $feature->seq_id;
				    my $ls   = abs( $s_start - $s_end );
				    return "$e ";
				  });
  }

  my $similar_length  = abs( $q_start - $q_end );
  my $s_simlar_length = abs( $s_start - $s_end );
  
  # Criterio de corte de HSB
  foreach my $currGroup (@group){
    if( $subject =~ $currGroup->{regexp}  ){
      if($e_value >  $currGroup->{e_threshold} ){
	next SUBJECT_ITERATION;
      } else {
	last;
      }
    }
  }

  # Identifying subject sequence group
  my $groupNum = getGroupNum( $subject );

  if( defined $coordTolerance && $groupNum != -1 ){
    # Coordenadas na ordem correta: start < end
    my ($q5Start, $q3End );

    if( $q_start < $q_end ){
      $q5Start = $q_start;
      $q3End   = $q_end;
    }else{
      $q5Start = $q_end;
      $q3End   = $q_start;
    }


    foreach my $currQryCov ( @{$qryCov[ $groupNum ]} ){

      # Do not paint if the region is already covered
      goto SUBJECT_ITERATION if( abs( $currQryCov->{start} - $q5Start ) < $coordTolerance  &&
				 abs( $currQryCov->{end}   - $q3End   ) < $coordTolerance ) ;    
    }

    # if this region was not covered till now
    my $ind = scalar( @{$qryCov[$groupNum]} );
    $qryCov[$groupNum][ $ind ]{start} = $q5Start;
    $qryCov[$groupNum][ $ind ]{end}   = $q3End;
  }

  my $feature = Bio::SeqFeature::Generic->new(-score     => $score,
					      #-seq_id      => "e=" . $e_value . " start:" . $s_start . " end:" . $s_end,
					      -seq_id      => "qc:$q_start-$q_end sc:$s_start-$s_end sl:$s_length s:$score l:$s_simlar_length e:$e_value %i:$perc_ident" ,
					      -display_name => $subject . " length:" . $s_length,
					      -start        => $q_start,
					      -end          => $q_end );

  $track->add_feature($feature);
}
#$progSubjIte->update( $lines );



#Saving picture

# Caso o usuario nao tenha fornecido o nome do arquivo
# de saida
if( ! defined( $outputFile ) ){
  $outputFile = $inputFile;
  if( !defined( $grpFile ) ){
    $outputFile .= "_$eThreForAll"
  }else{
    foreach my $currGroup (@group){
	$outputFile .= "_" . $currGroup->{name} . "_" . $currGroup->{color} .
		      "_" . $currGroup->{e_threshold}
    }
  }
  $outputFile .= "_c_" . $coordTolerance if( defined( $coordTolerance ) );
  $outputFile .= ".svg";
}
open OUTPUT, ">" . $outputFile;
print OUTPUT $panel->svg;
close OUTPUT;

print "\nDisplaying graph...\n";
#system ("gimp $outputFile");

exit();

sub getGroupNum{
  my ($seqId) = @_;

  for ( my $indGroup = 0;
	$indGroup < scalar( @group );
	$indGroup++ ) {
    return $indGroup 
      if ( $seqId =~ /$group[$indGroup]->{regexp}/ );
  }
  return -1;
}
