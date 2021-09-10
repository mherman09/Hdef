#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Text::CSV;


# Usage statement
sub usage {
    print "Usage: parse_usgs_query.pl --input_file FILE [...options...]\n";
    print "-i|--input_file FILE            Name of input CSV file\n";
    print "-o|--origin_time                Origin time\n";
    print "-l|--location                   Location\n";
    print "-m|--magnitude                  Magnitude\n";
    print "-t|--tensor                     Moment tensor\n";
    print "-p|--priority MT_TYP            Select priority MT type (Mww,Mwc,Mwr,Mwb,duputel_Mww)\n";
    die;
}

# Parse the command line
my $input_file = '';
my $iWantOriginTime = 0;
my $iWantLocation = 0;
my $iWantMagnitude = 0;
my $iWantMomentTensor = 0;
my $priorityMTType = 'Mww';
GetOptions(
    'input_file=s' => \$input_file,
    'origin_time' => \$iWantOriginTime,
    'location' => \$iWantLocation,
    'magnitude' => \$iWantMagnitude,
    'tensor' => \$iWantMomentTensor,
    'priority=s' => \$priorityMTType
);


# Open the input file
open my $fh, "<", $input_file or &usage;


# Parse the input file
my $csv = Text::CSV->new ({
    binary    => 1, # Allow special character. Always set this
    auto_diag => 1, # Report irregularities immediately
});


# Initialize variables
my $i_time = -1;
my $i_longitude = -1;
my $i_latitude = -1;
my $i_depth = -1;
my $i_magnitude = -1;
my $i_us_Mww = -1;
my $i_us_Mwb = -1;
my $i_us_Mwr = -1;
my $i_us_Mwc = -1;
my $i_duputel_Mww = -1;

# First line of the file is a header with field information
# Get the indices of all variables
my $row = $csv->getline($fh);
my $i = 0;
for (@$row) {
    if ($_ =~ "longitude") {
        $i_longitude = $i;
        #print "i_longitude=$i_longitude\n";
    } elsif ($_ =~ "latitude") {
        $i_latitude = $i;
        #print "i_latitude=$i_latitude\n";
    } elsif ($_ =~ "depth") {
        $i_depth = $i;
        #print "i_depth=$i_depth\n";
    } elsif ($_ =~ "magnitude") {
        $i_magnitude = $i;
        #print "i_magnitude=$i_magnitude\n";
    } elsif ($_ =~ "time") {
        $i_time = $i;
        #print "i_time=$i_time\n";
    } elsif ($_ =~ "us_Mww_mrr") {
        $i_us_Mww = $i;
        #print "i_us_Mww=$i_us_Mww\n";
    } elsif ($_ =~ "us_Mwr_mrr") {
        $i_us_Mwr = $i;
        #print "i_us_Mwr=$i_us_Mwr\n";
    } elsif ($_ =~ "us_Mwb_mrr") {
        $i_us_Mwb = $i;
        #print "i_us_Mwb=$i_us_Mwb\n";
    } elsif ($_ =~ "us_Mwc_mrr") {
        $i_us_Mwc = $i;
        #print "i_us_Mwc=$i_us_Mwc\n";
    } elsif ($_ =~ "duputel_Mww_mrr") {
        $i_duputel_Mww = $i;
        #print "i_us_Mwb=$i_duputel_Mww\n";
    }
    $i = $i + 1;
}


# Go through the rest of the file and print the requested quantities
my $iWantSomething = 0;
while (my $row = $csv->getline ($fh)) {
    my $ot = $row->[$i_time];
    my $lon = sprintf "%10.4f", $row->[$i_longitude];
    my $lat = sprintf "%8.4f", $row->[$i_latitude];
    my $dep = sprintf "%6.2f", $row->[$i_depth];
    my $mag = sprintf "%6.2f", $row->[$i_magnitude];
    my $i_mt = -1;
    my $mt_type = "";
    my @mij = (0,0,0,0,0,0);

    # Build the output string
    my $output = "";
    if ($iWantOriginTime) {
        if ($i_time >= 0) {
            $output = $output."$ot ";
        } else {
            print "parse_usgs_query.pl: requested origin time but \"time\" was not found in the header\n";
            die;
        }
        $iWantSomething = 1;
    }
    if ($iWantLocation) {
        if ($i_longitude >= 0) {
            $output = $output."$lon ";
        } else {
            print "parse_usgs_query.pl: requested location but \"longitude\" was not found in the header\n";
            die;
        }
        if ($i_latitude >= 0) {
            $output = $output."$lat ";
        } else {
            print "parse_usgs_query.pl: requested location but \"latitude\" was not found in the header\n";
            die;
        }
        if ($i_depth >= 0) {
            $output = $output."$dep ";
        } else {
            print "parse_usgs_query.pl: requested location but \"depth\" was not found in the header\n";
            die;
        }
        $iWantSomething = 1;
    }
    if ($iWantMagnitude) {
        if ($i_magnitude >= 0) {
            $output = $output."$mag ";
        } else {
            print "parse_usgs_query.pl: requested magnitude but \"magnitude\" was not found in the header\n";
            die;
        }
        $iWantSomething = 1;
    }
    if ($iWantMomentTensor) {
        # Specifically defined moment tensor priority
        if ($priorityMTType ne "") {
            if ($priorityMTType eq "Mww" && $i_us_Mww >= 0 && $row->[$i_us_Mww] ne "") {
                $i_mt = $i_us_Mww;
                $mt_type = "Mww";
            } elsif ($priorityMTType eq "Mwc" && $i_us_Mwc >= 0 && $row->[$i_us_Mwc] ne "") {
                $i_mt = $i_us_Mwc;
                $mt_type = "Mwc";
            } elsif ($priorityMTType eq "Mwr" && $i_us_Mwr >= 0 && $row->[$i_us_Mwr] ne "") {
                $i_mt = $i_us_Mwr;
                $mt_type = "Mwr";
            } elsif ($priorityMTType eq "Mwb" && $i_us_Mwb >= 0 && $row->[$i_us_Mwb] ne "") {
                $i_mt = $i_us_Mwb;
                $mt_type = "Mwb";
            } elsif ($priorityMTType eq "duputel_Mww" && $i_duputel_Mww >= 0 && $row->[$i_duputel_Mww] ne "") {
                $i_mt = $i_duputel_Mww;
                $mt_type = "duputel_Mww";
            }
        }

        # Default moment tensor priority list
        if ($mt_type eq "") {
            if ($i_us_Mww >= 0 && $row->[$i_us_Mww] ne "") {
                $i_mt = $i_us_Mww;
                $mt_type = "Mww";
            } elsif ($i_us_Mwc >= 0 && $row->[$i_us_Mwc] ne "") {
                $i_mt = $i_us_Mwc;
                $mt_type = "Mwc";
            } elsif ($i_us_Mwr >= 0 && $row->[$i_us_Mwr] ne "") {
                $i_mt = $i_us_Mwr;
                $mt_type = "Mwr";
            } elsif ($i_us_Mwb >= 0 && $row->[$i_us_Mwb] ne "") {
                $i_mt = $i_us_Mwb;
                $mt_type = "Mwb";
            } elsif ($i_duputel_Mww >= 0 && $row->[$i_duputel_Mww] ne "") {
                $i_mt = $i_duputel_Mww;
                $mt_type = "duputel_Mww";
            } else {
                $output = $output."no_MT ";
            }
        }
        if ($i_mt >= 0) {
            $mij[0] = sprintf "%12.4e", $row->[$i_mt];
            $mij[1] = sprintf "%12.4e", $row->[$i_mt+1];
            $mij[2] = sprintf "%12.4e", $row->[$i_mt+2];
            $mij[3] = sprintf "%12.4e", $row->[$i_mt+3];
            $mij[4] = sprintf "%12.4e", $row->[$i_mt+4];
            $mij[5] = sprintf "%12.4e", $row->[$i_mt+5];
            $output = $output."@mij $mt_type ";
        }
        $iWantSomething = 1;
    }
    if ($iWantSomething eq 1) {
        print "$output\n";
    } else {
        print "parse_usgs_query.pl: no output defined\n";
        &usage;
    }
}


# Close the file
close $fh;
