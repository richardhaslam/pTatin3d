#!/usr/bin/perl

 use Cwd;

#my $co_name = '/Users/dmay/codes/ptat3d-rh-riftspm/ptat3d-subduction-dbg';
my $co_name = getcwd . '/..';

`svn info $co_name/src > .SVN_INFO`;
`svn log $co_name/src --limit 1 > .SVN_LOG`;

$file = '.SVN_INFO';		# Name the file
open(INFO, $file);		# Open the file
@lines_info = <INFO>;		# Read it into an array
close(INFO);			# Close the file

$file = '.SVN_LOG';		# Name the file
open(INFO, $file);		# Open the file
@lines_log = <INFO>;		# Read it into an array
close(INFO);			# Close the file

unlink '.SVN_INFO';
unlink '.SVN_LOG';

chomp( @lines_info[3] );
chomp( @lines_info[4] );
chomp( @lines_log[1]  );

open (MYFILE, '>ptatin_svn_info.h');

print MYFILE "\n";
print MYFILE "#ifndef __ptatin_svn_info_h__\n";
print MYFILE "#define __ptatin_svn_info_h__\n";
print MYFILE "\n";
print MYFILE "  const char PTATIN_SVN_REPO_UUID[] = \"@lines_info[3]\"\; \n";
print MYFILE "  const char PTATIN_SVN_REPO_REV[] = \"@lines_info[4]\"\; \n";
print MYFILE "  const char PTATIN_SVN_REPO_LAST_CHANGE[] = \"Log: @lines_log[1]\"\; \n";
print MYFILE "\n";
print MYFILE "#endif \n";
print MYFILE "\n";

close (MYFILE);