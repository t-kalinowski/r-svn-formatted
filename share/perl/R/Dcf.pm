package R::Dcf;

=head1 NAME
    
  R::Dcf - parse files in dcf format
    
=head1 SYNOPSIS

  use R::Dcf;
  $rdcf = R::Dcf->new( "foo" );

  @keys = keys( %{ $rdcf } );     # returns an array of the fields in the file
  $field1 = $rdcf->{ $keys[0] };  # returns the first field

=cut
          
use strict;
use Carp;
use FileHandle;


sub new {

    my ($class, $file) = @_;
    $class = ref($class) || $class;

    my $fh = new FileHandle "< $file" or croak "open($file): $!\n";
    my $self = {};
    my $field = "";
    while(<$fh>){
	s/\r/ /g;
	s/\s*$//g;
	if(/^\s+/){
	    if($field){
		s/^\s+/\n/;
		$self->{ $field } .= $_;
	    }
	    else{
		die "Malformed DCF file (line $.)\n";
	    }
	}
	else{
	    if(/^(\S+):\s*(.*\S)\s*$/){
		$field=$1;
		$self->{ $field } .= $2;
	    }
	    else{
		die "Malformed DCF file (line $.)\n";
	    }
	}
    }	
    $fh->close;
    bless $self, $class;
}

1;




