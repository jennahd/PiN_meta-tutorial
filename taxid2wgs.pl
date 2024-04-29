#!/usr/bin/env perl
# $Id: taxid2wgs.pl,v 1.6 2017/06/23 13:53:09 camacho Exp $
use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Long;
use Pod::Usage;
use autodie;

use constant URL => "https://www.ncbi.nlm.nih.gov/blast/BDB2EZ/taxid2wgs.cgi";
use constant INCLUDE => "INCLUDE_TAXIDS";
use constant EXCLUDE => "EXCLUDE_TAXIDS";
use constant USER_AGENT => "taxid2wgs.pl/1.1";
use constant EXTN => ".nvl";

my @includes = ();
my @excludes = ();
my $url_api_ready = 0;
my ($title, $alias_file);
my $verbose = 0;
my $help_requested = 0;
GetOptions("excludes=i"     => \@excludes,
           "title=s"        => \$title,
           "alias_file=s"   => \$alias_file,
           "verbose|v+"     => \$verbose,
           "url_api_ready"  => \$url_api_ready,
           "help|?"         => \$help_requested) || pod2usage(2);
pod2usage(-verbose=>2) if ($help_requested);
pod2usage(-msg => "Must provide taxids", -verbose=>2) unless (@ARGV);
pod2usage(-msg => "Alias file cannot be empty", -verbose=>2) 
    if (defined($alias_file) and (length($alias_file) == 0));
pod2usage(-msg => "Title must be provided with alias_file", -verbose=>2) 
    if (defined($title) and not defined($alias_file));
pod2usage(-msg => "Title cannot be empty", -verbose=>2) 
    if (defined($title) and (length($title) == 0));

my @projects = &get_wgs_projects();
die "Failed to get taxids\n" unless (@projects);

if (defined $alias_file) {
    map { s,WGS_VDB://,, } @projects;
    &print_alias($alias_file, $title, \@projects);
} else {
    map { s,WGS_VDB://,, } @projects unless ($url_api_ready);
    print join("\n", @projects);
}

sub get_wgs_projects
{
    my $url = URL . "?" . INCLUDE . "=" . join(",", @ARGV);
    $url .= "&" . EXCLUDE . "=" . join(",", @excludes) if scalar (@excludes);
    my $ua = LWP::UserAgent->new;
    $ua->agent(USER_AGENT);
    $ua->show_progress(1) if ($verbose >= 3);
    if ($verbose >= 2) {
        $ua->add_handler("request_send", sub { shift->dump; return; });
        $ua->add_handler("response_done", sub { shift->dump; return; });
    }
    my $req = HTTP::Request->new(POST => $url);
    my $res = $ua->request($req);
    print $res->decoded_content if $verbose;
    my @retval;
    if ($res->is_success) {
        @retval = split(/\s+/, $res->content);
    } else {
        print STDERR $res->status_line, "\n";
    }
    return @retval;
}

sub print_alias
{
    use Time::localtime;
    my $basename = shift;
    my $title = shift;
    my $projects = shift;
    open(my $out, ">", $basename . EXTN);
    print $out "#\n# Alias file created by " . USER_AGENT;
    print $out " on " . ctime(time) . "\n";
    print $out "TITLE $title\n" if (defined $title and length $title);
    print $out "VDBLIST " . join(" ", @$projects) . "\n";
    close ($out);
}

__END__

=head1 NAME

B<taxid2wgs.pl> - Retrieve WGS projects for given NCBI taxonomy IDs

=head1 SYNOPSIS

taxid2wgs.pl [options] <taxid1> [ <taxid2> ... <taxidN> ]

=head1 ARGUMENTS

=over

=item B<-exclude>

Accepts taxids to exclude from the output (default: None)

=item B<-alias_file>

File base name (no extension) to save results into an alias file.

=item B<-title>

Title to include in the generated alias file. Required if alias_file is provided.

=item B<-url_api_ready>

Produce output that can be used in the NCBI URL API (default: false)

=item B<-verbose>, B<-v>

Produce verbose output, can be specified multiple times for increased verbosity (default: false)

=item B<-help>, B<-?>

Displays this man page.

=back

=head1 AUTHOR

Christiam Camacho (camacho@ncbi.nlm.nih.gov)

=cut

