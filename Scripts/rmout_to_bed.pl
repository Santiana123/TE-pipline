#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($input_file, $bed_file);
GetOptions(
    'input=s' => \$input_file,
    'bed=s'    => \$bed_file,
) or die "Error in command line arguments";

die "Usage: perl extract_bed.pl -input <input_file> -bed <bed_file>\n" unless $input_file && $bed_file;

open(my $in_fh, '<', $input_file) or die "Could not open file '$input_file' $!";
open(my $bed_fh, '>', $bed_file) or die "Could not open file '$bed_file' $!";

while (<$in_fh>) {
    chomp;
    my @data = split(/\t/);  # 使用制表符分割
    next if (@data < 7);  # 确保有足够的列

    next unless grep { /Unknown/ } @data;

    my ($chr, $start, $end) = ($data[4], $data[5], $data[6]);
    
    $start -= 1;  # 转换为0-based

    print $bed_fh join("\t", $chr, $start, $end, '.', 0, '+') . "\n";
}

close($in_fh);
close($bed_fh);

print "BED file has been created: '$bed_file'.\n";

