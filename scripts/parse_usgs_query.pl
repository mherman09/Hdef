#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Text::CSV;


# Usage statement
sub usage {
    print STDERR "Usage: parse_usgs_query.pl --input_file FILE [...options...]\n";
    print STDERR "-i|--input_file FILE            Name of input CSV file\n";
    print STDERR "-o|--origin_time                Origin time\n";
    print STDERR "-l|--location                   Location\n";
    print STDERR "-m|--magnitude                  Magnitude\n";
    print STDERR "-t|--tensor                     Moment tensor\n";
    print STDERR "-f|--focal_mechanism            Focal mechanism\n";
    print STDERR "-p|--priority MT_TYP            Select priority MT type (Mww,Mwc,Mwr,Mwb,duputel_Mww)\n";
    die;
}

# Parse the command line
my $input_file = '';
my $iWantOriginTime = 0;
my $iWantLocation = 0;
my $iWantMagnitude = 0;
my $iWantMomentTensor = 0;
my $priorityMTType = 'Mww';
my $iWantFocalMechanism = 0;
GetOptions(
    'input_file=s' => \$input_file,
    'origin_time' => \$iWantOriginTime,
    'location' => \$iWantLocation,
    'magnitude' => \$iWantMagnitude,
    'tensor' => \$iWantMomentTensor,
    'focal_mechanism' => \$iWantFocalMechanism,
    'priority=s' => \$priorityMTType
);


# Check for the input file
if (! -f $input_file) {
    print STDERR "parse_usgs_query.pl: could not find seismicity file named \"$input_file\"\n";
    &usage;
}

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
my $i_us_Mww = -1; # Moment tensors
my $i_us_Mwb = -1;
my $i_us_Mwr = -1;
my $i_us_Mwc = -1;
my $i_duputel_Mww = -1;
my $i_nc_TMTS = -1;
my $i_nc_TDMT = -1;
my $i_nc_Mw = -1;
my $i_nn_Mw = -1;
my $i_nc_Mww = -1;
my $i_gcmt_np = -1; # Focal mechanisms
my $i_us_np = -1;
my $i_nc_np = -1;

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
    } elsif ($_ =~ "magnitude" || $_ eq "mag") {
        $i_magnitude = $i;
        #print "i_magnitude=$i_magnitude\n";
    } elsif ($_ =~ "time") {
        $i_time = $i;
        #print "i_time=$i_time\n";         # Moment tensors
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
    } elsif ($_ =~ "nc_TMTS_mrr") {
        $i_nc_TMTS = $i;
    } elsif ($_ =~ "nc_TDMT_mrr") {
        $i_nc_TDMT = $i;
    } elsif ($_ =~ "nn_Mw_mrr") {
        $i_nn_Mw = $i;
    } elsif ($_ =~ "nc_Mww_mrr") {
        $i_nc_Mww = $i;
    } elsif ($_ =~ "nc_Mw_mrr") {
        $i_nc_Mw = $i;
    } elsif ($_ =~ "gcmt_np1_strike") {    # Focal mechanisms
        $i_gcmt_np = $i;
    } elsif ($_ =~ "us_np1_strike") {
        $i_us_np = $i;
    } elsif ($_ =~ "nc_np1_strike") {
        $i_nc_np = $i;
    }
    $i = $i + 1;
}


# Go through the rest of the file and print the requested quantities
$i = 2;
my $iWantSomething = 0;
while (my $row = $csv->getline ($fh)) {
    my $ot = $row->[$i_time];
    my $lon = sprintf "%10.4f", $row->[$i_longitude];
    my $lat = sprintf "%8.4f", $row->[$i_latitude];
    my $dep = 0.0;
    if ($row->[$i_depth] ne "") {
        $dep = sprintf "%6.2f", $row->[$i_depth];
    } else {
        print STDERR "parse_usgs_query.pl: [WARNING] line $i has blank depth\n";
    }
    my $mag = sprintf "%6.2f", $row->[$i_magnitude];
    my $i_mt = -1;
    my $mt_type = "";
    my @mij = (0,0,0,0,0,0);
    my $i_fm = -1;
    my $fm_type = "";
    my @fm = (0,0,0,0,0,0);

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
            } elsif ($i_nc_Mw >= 0 && $row->[$i_nc_Mw] ne "") {
                $i_mt = $i_nc_Mw;
                $mt_type = "nc_Mw";
            } elsif ($i_nc_Mww >= 0 && $row->[$i_nc_Mww] ne "") {
                $i_mt = $i_nc_Mww;
                $mt_type = "nc_Mww";
            } elsif ($i_nn_Mw >= 0 && $row->[$i_nn_Mw] ne "") {
                $i_mt = $i_nn_Mw;
                $mt_type = "nn_Mw";
            } elsif ($i_nc_TDMT >= 0 && $row->[$i_nc_TDMT] ne "") {
                $i_mt = $i_nc_TDMT;
                $mt_type = "nc_TDMT";
            } elsif ($i_nc_TMTS >= 0 && $row->[$i_nc_TMTS] ne "") {
                $i_mt = $i_nc_TMTS;
                $mt_type = "nc_TMTS";
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
    if ($iWantFocalMechanism) {
        # Default focal mechanism priority list
        if ($mt_type eq "") {
            if ($i_us_np >= 0 && $row->[$i_us_np] ne "") {
                $i_fm = $i_us_np;
                $fm_type = "us";
            } elsif ($i_gcmt_np >= 0 && $row->[$i_gcmt_np] ne "") {
                $i_fm = $i_gcmt_np;
                $fm_type = "gcmt";
            } elsif ($i_nc_np >= 0 && $row->[$i_nc_np] ne "") {
                $i_fm = $i_nc_np;
                $fm_type = "nc";
            } else {
                $output = $output."no_FM ";
            }
        }
        if ($i_fm >= 0) {
            $fm[0] = sprintf "%12.4e", $row->[$i_fm];
            $fm[1] = sprintf "%12.4e", $row->[$i_fm+1];
            $fm[2] = sprintf "%12.4e", $row->[$i_fm+2];
            $fm[3] = sprintf "%12.4e", $row->[$i_fm+3];
            $fm[4] = sprintf "%12.4e", $row->[$i_fm+4];
            $fm[5] = sprintf "%12.4e", $row->[$i_fm+5];
            $output = $output."@fm $fm_type ";
        }
        $iWantSomething = 1;
    }
    if ($iWantSomething eq 1) {
        print "$output\n";
    } else {
        print STDERR "parse_usgs_query.pl: no output defined\n";
        &usage;
    }
    $i = $i + 1;
}


# Close the file
close $fh;
