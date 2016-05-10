#!/usr/bin/perl

# Adds a copyright header to all .c and .h files.
# The header block is placed between the comments /*@ @*/, this allows us
# replace existing blocks with another variation

## NOTE: After running this script you will have to possibly update the svn propset keyword Id
## To do this, do the following
##    for i in `find . -name "*.c"`; do svn propset svn:keywords "Id" $i; done
##

sub wanted_c_h_files;
sub process_file;
sub print_header;
sub print_header_2008;
sub trim($);


use File::Find;

# 1) Generate a list of all *.c and *.h files (full path to file) beneath current directory

## Set equal to 'true' when debugging. Thus will make a copy of each file with .bup extension
## To remove all the buckups do > for i in `find . -name "*.bup"`; do rm -f $i; done
$paranoid = 'false';

find( \&wanted_c_h_files, './' );
#print @file_list;

sub wanted_c_h_files 
{
    $full_path = $File::Find::dir;
    $full_path_to_file = $File::Find::name;
    $cur_file = $_;
  
    # check extension of $cur_file, if .c or .h, then add it
    if( $cur_file =~ m/\.c$/ || $cur_file =~ m/\.h$/) {
    #    print "  Found (.c , .h) $full_path_to_file \n";
        process_file( $full_path_to_file, $cur_file );
    }
}

#
#  We assume there is no data on the same line following the close tag
#
sub process_file
{
    my $file     = $_[0];
    my $filename = $_[1];

    my @copy_file;
    my $n = 0;
    my $state = 0;
    
    
    # scan file and ONLY copy lines which are NOT between the special tags
    open( DAT, $filename ) || die( "Could not open file ($filename) ! \n" );
    @raw_data=<DAT>;
    foreach $line (@raw_data) {
        $tmp_line = $line;
        # look for open tag at start of line
        if( $line =~ m/^\/\*\@/ ) {
         #   print "Found an open tag \n";
            $state++;
        }

        # Could check for close at end of line in one hit....
        if( $line =~ m/\@\*\/$/ ) {
        #    print "Verified close tag is at the end of the line \n";
            $state--;
             next;
        }
 
        if( $state == 0 ) {
            $copy_file[ $n ] = $tmp_line;
            $n++;
        }
    }
    close( DAT );

    if( $state > 0 ) {
        print "===============================================================================\n";
        print "*************                       WARNING                       *************\n";
        print "  It appears that the open tag was never matched by a close tag, \n";
        print "  thus file content may be missing. This indicates an error.\n";
        print "  Possible causes i)  the close tag @*/ was not used \n";
        print "                  ii) close tag did not appear at the END of the line. \n";
        print "  AddReplaceCopyrightHeader was NOT applied to file $filename \n";
        print "===============================================================================\n\n";
        
        return;
    }
    if( $state < 0 ) {
        print "===============================================================================\n";
        print "*************                       WARNING                       *************\n";
        print "  It appears that the close tag was never matched by an open tag, \n";
        print "  thus file content may be missing. This indicates an error.\n";
        print "  Possible causes i)  the open tag /*@ was not used \n";
        print "                  ii) open tag did not appear at the BEGINNING of the line. \n";
        print "  AddReplaceCopyrightHeader was NOT applied to file $filename \n";
        print "===============================================================================\n\n";
        
        return;
    }
    
    # Make a backup (paranoid)
    if( $paranoid eq 'true' ) {
        $bup = $filename . '.bup';
        `cp $filename $bup`;
    }

    # create a new file with the same name 
    # Add header at the of the file

		print_header_2012( 'tmp.file', 'pTatin3d', $filename );
	

    # Dump everything else
    open( DAT, ">>tmp.file" ) || die( "Could not open file ($filename) ! \n" );
    foreach $line (@copy_file) {
        print DAT "$line";
    }
    close( DAT );

#    print "    Moving tmp.file -> $filename \n";
    `mv tmp.file $filename`;

    return ;    
}

sub print_header_2012
{
    my $filename = $_[0];
    my $proj = $_[1];
    my $name = $_[2];
	
    my $sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst;
	
    ## Get the time-date information
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $year += 1900;
	
    open( DAT, ">$filename" ) || die( "Could not open file ($filename) ! \n" );
	
    print DAT "/*@";
    
    print DAT " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may\@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    $proj
 **    filename:   $name
 **
 **
 **    $proj is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    $proj is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with $proj. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
		
    print DAT " @*/\n";
	
    close( DAT );
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

