#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($input_file, $output_file);
GetOptions(
    'input=s'  => \$input_file,
    'output=s' => \$output_file,
) or die "Error in command line arguments";

die "Usage: perl all-types-stat.pl -input <input_file> -output <output_file>\n" unless $input_file && $output_file;

my %type_count;
my %type_length;

open(my $in_fh, '<', $input_file) or die "Could not open file '$input_file' $!";
while (<$in_fh>) {
    chomp;
    my @data = split(/\s+/);
    next if (@data < 11);  # 确保有足够的列

    my $type = $data[10];
    my $start = $data[5];
    my $end = $data[6];

    $type_count{$type}++;
    $type_length{$type} += ($end - $start);
}
close($in_fh);

open(my $out_fh, '>', $output_file) or die "Could not open file '$output_file' $!";
foreach my $type (sort keys %type_count) {
    print $out_fh "$type\t$type_count{$type}\t$type_length{$type}\n";
}
close($out_fh);

print "Statistics have been written to '$output_file'.\n";
