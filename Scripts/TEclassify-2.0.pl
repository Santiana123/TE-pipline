#!/usr/bin/perl 

use Getopt::Std;
getopts "i:g:u:";


if ((!defined $opt_i)|| (!defined $opt_g)||(!defined $opt_u) ) {
    die "************************************************************************
    Usage: perl $0 -i repeatmasker.out -u unknown.TE.fasta.lib -g genome_size
           -i : repeatmasker.out
           -g : genome_size, could be [k|K|m|M|g|G]
           -u : unknown.TE.fasta.lib, results from TEclass
           -o : result.out
************************************************************************\n";
}

my $g_size = $opt_g;
if($g_size=~/[k|K]/){
	$g_size =~ s/[k|K]//g;
	$g_size = $g_size * 1000;
}elsif($g_size=~/[m|M]/){
	$g_size =~ s/[m|M]//g;
	$g_size = $g_size * 1000000;	
}elsif($g_size=~/g|G/){		
	$g_size =~ s/[g|G]//g;
	$g_size = $g_size * 1000000000;	
	}

my %infordb;
my %alldb;
open(IN, $opt_i) or die"";
<IN>;
while(<IN>){
	chomp;
	next if(/score/);
	my @data = split(/\s+/,$_);
	my $num  = @data;
	next if($num < 5);
	my $family = $data[-5] if($data[-1] ne "*");
	   $family = $data[-6] if($data[-1] eq "*");
	my $l    = abs($data[-9] - $data[-10]) + 1 if($data[-1] ne "*");
	   $l    = abs($data[-10] - $data[-11]) + 1 if($data[-1] eq "*");
	next if($family eq " ");
	next if($family eq "class/family");
	next if(!defined $family);
	$infordb{$family}->{'size'} += $l;
  $infordb{$family}->{'num'}  += 1;  
  my $contig = ($data[-1] ne "*")?$data[-11]:$data[-12];
  my $a      = ($data[-1] ne "*")?$data[-10]:$data[-11];
  my $b      = ($data[-1] ne "*")?$data[-9]:$data[-10];
	}
close IN;

my %tRepdb;
open(INN, "tandem.all.bed") or die"";
while(<INN>){
	chomp;
	my @data = split(/\s+/,$_);
	$infordb{'tRep'}->{'num'}++;
	foreach my $i($data[1]..$data[2]){
		$alldb{$contig}->{$i}++;
		$tRepdb{$contig}->{$i}++;
		}
	}
close INN;

foreach my $contig(keys %tRepdb){
	$infordb{'tRep'}->{'size'} += keys %{$tRepdb{$contig}}; 
	}



open(OUT, "> tmp.TE.out") or die"";
foreach my $family (sort keys %infordb){
	print OUT "$family	$infordb{$family}->{'num'}	$infordb{$family}->{'size'}\n";
	}
close OUT;

open(IN, "tmp.TE.out") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	if(/LTR/){
		$infordb{'LTR'}->{'num'}  += $data[1];
		$infordb{'LTR'}->{'size'} += $data[2];
	}elsif(/DNA/){
		$infordb{'DNA'}->{'num'}  += $data[1];
		$infordb{'DNA'}->{'size'} += $data[2];		
		}
	
	if(/SINE/){
		$infordb{'SINE'}->{'num'}  += $data[1];
		$infordb{'SINE'}->{'size'} += $data[2];
	}elsif(/LINE/){
		$infordb{'LINE'}->{'num'}  += $data[1];
		$infordb{'LINE'}->{'size'} += $data[2];		
	}elsif($_=~/rRNA/ or $_=~/snRNA/){
		$infordb{'sRNA'}->{'num'}  += $data[1];
		$infordb{'sRNA'}->{'size'} += $data[2];			
	}elsif($_=~/Simple_repeat/ or $_=~/Low complexity/){
		$infordb{'tRep'}->{'num'}  += $data[1];
		$infordb{'tRep'}->{'size'} += $data[2];		
	}elsif(/Unknown/){
		$infordb{'unKnw'}->{'num'}  += $data[1];
		$infordb{'unKnw'}->{'size'} += $data[2];		
    next;	
	}elsif(/Copia/){
		$infordb{'copia'}->{'num'}  += $data[1];
		$infordb{'copia'}->{'size'} += $data[2];				
	}elsif(/Gypsy/){
		$infordb{'gypsy'}->{'num'}  += $data[1];
		$infordb{'gypsy'}->{'size'} += $data[2];
        }elsif(/DNA\/DTT/ or /MITE\/DTT/ or /DNA\/TcMar/ or /DNA\/TcMar-Pogo/){
                $infordb{'TcMar'}->{'num'}  += $data[1];
                $infordb{'TcMar'}->{'size'} += $data[2];
        }elsif(/DNA\/DTC/ or /MITE\/DTC/ or /DNA\/En-Spm/){
                $infordb{'EnSpm/CACTA'}->{'num'}  += $data[1];
                $infordb{'EnSpm/CACTA'}->{'size'} += $data[2];
        }elsif(/DNA\/Harbinger/ or /DNA\/DTH/ or /MITE\/DTH/ or /PIF-Harbinger/){
                $infordb{'PIF-Harbinger'}->{'num'}  += $data[1];
                $infordb{'PIF-Harbinger'}->{'size'} += $data[2];
        }elsif(/DNA\/DTM/ or /MITE\/DTM/ or /DNA\/MuDR/ or $_=~/Mu/ or $_=~/MU/){
                $infordb{'MuDR/Mutator'}->{'num'}  += $data[1];
                $infordb{'MuDR/Mutator'}->{'size'} += $data[2];
        }elsif(/DNA\/Helitron/ or /RC\/Helitron/){
                $infordb{'Helitron'}->{'num'}  += $data[1];
                $infordb{'Helitron'}->{'size'} += $data[2];
        }elsif(/hAT/ or /DNA\/hAT/ or /DNA\/hAT-Ac/ or /DNA\/DTA/ or /MITE\/DTA/){
                $infordb{'hAT'}->{'num'}  += $data[1];
                $infordb{'hAT'}->{'size'} += $data[2];
        }elsif(/Retroposon/){
                $infordb{'otherRetro'}->{'num'}  += $data[1];
                $infordb{'otherRetro'}->{'size'} += $data[2];
		}
	}
close IN;

open(IN, $opt_u) or die"";
<IN>;
$/='>';
while(<IN>){
	chomp;
	my ($TE,$seq) = split(/\n/,$_,2);
	my $len_TE    = length $seq;
	if($TE =~ /result:\s+DNA/){
		$infordb{'DNA'}->{'num'}  += 1;
		$infordb{'DNA'}->{'size'} += $len_TE;			
	}elsif($TE =~ /result:\s+SINE/){
		$infordb{'SINE'}->{'num'}  += 1;
		$infordb{'SINE'}->{'size'} += $len_TE;		
	}elsif($TE =~ /result:\s+LINE/){
		$infordb{'LINE'}->{'num'}  += 1;
		$infordb{'LINE'}->{'size'} += $len_TE;			
	}elsif($TE =~ /result:\s+Retro/ or $TE =~ /result:\s+nonLTR/){
		$infordb{'otherRetro'}->{'num'}  += 1;
		$infordb{'otherRetro'}->{'size'} += $len_TE;			
	}elsif($TE =~ /result:\s+LTR/){
		$infordb{'LTR'}->{'num'}  += 1;
		$infordb{'LTR'}->{'size'} += $len_TE;		
	}elsif($TE =~ /result:\s+unclear/){
		$infordb{'unKnw'}->{'num'}  += 1;
		$infordb{'unKnw'}->{'size'} += $len_TE;			
		}
	}
close IN;

my $non_LTR_num  = $infordb{'SINE'}->{'num'} + $infordb{'LINE'}->{'num'} + $infordb{'otherRetro'}->{'num'};
my $non_LTR_size = $infordb{'SINE'}->{'size'} + $infordb{'LINE'}->{'size'} + $infordb{'otherRetro'}->{'size'};

my $retro_num    = $infordb{'LTR'}->{'num'} + $non_LTR_num;
my $retro_size   = $infordb{'LTR'}->{'size'} + $non_LTR_size;

$total_repeat = $retro_size + $infordb{'DNA'}->{'size'} + $infordb{'Helitron'}->{'size'} + $infordb{'tRep'}->{'size'} + $infordb{'unKnw'}->{'size'};
$count = $retro_num + $infordb{'DNA'}->{'num'} + $infordb{'Helitron'}->{'num'} + $infordb{'tRep'}->{'num'} +  $infordb{'unKnw'}->{'num'};

$percent  =  sprintf("%.2f",$total_repeat/$g_size * 100);
$r_per    = 100.00;
print "Total repeat fraction    $count  $total_repeat   $r_per  $percent\n";

$percent    = sprintf("%.2f",$retro_size/$g_size * 100);
$r_per      = sprintf("%.2f",$retro_size/$total_repeat * 100);
print "Class I: Retroelement	$retro_num	$retro_size	$r_per	$percent\n";
$percent    = $infordb{'LTR'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'LTR'}->{'size'}/$total_repeat * 100);
print "	LTR Retrotransposon	$infordb{'LTR'}->{'num'}	$infordb{'LTR'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'copia'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'copia'}->{'size'}/$total_repeat * 100);
print "		Ty1/Copia	$infordb{'copia'}->{'num'}	$infordb{'copia'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'gypsy'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'gypsy'}->{'size'}/$total_repeat * 100);
print "		Ty3/Gypsy	$infordb{'gypsy'}->{'num'}	$infordb{'gypsy'}->{'size'}	$r_per	$percent\n";

my $other_LTR_num   = $infordb{'LTR'}->{'num'} - $infordb{'copia'}->{'num'} - $infordb{'gypsy'}->{'num'} ;
my $other_LTR_size  = $infordb{'LTR'}->{'size'} - $infordb{'copia'}->{'size'} - $infordb{'gypsy'}->{'size'};
$percent    = $other_LTR_size/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$other_LTR_size/$total_repeat * 100);
print "		Other	$other_LTR_num	$other_LTR_size	$r_per	$percent\n";
$percent    = $non_LTR_size/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$non_LTR_size/$total_repeat * 100);
print "	non-LTR Retrotransposon	$non_LTR_num	$non_LTR_size	$r_per	$percent\n";
$percent    = $infordb{'LINE'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'LINE'}->{'size'}/$total_repeat * 100);
print "		LINE	$infordb{'LINE'}->{'num'}	$infordb{'LINE'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'SINE'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'SINE'}->{'size'}/$total_repeat * 100);
print "		SINE	$infordb{'SINE'}->{'num'}	$infordb{'SINE'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'otherRetro'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'otherRetro'}->{'size'}/$total_repeat * 100);
print "	Other	$infordb{'otherRetro'}->{'num'}	$infordb{'otherRetro'}->{'size'}	$r_per	$percent\n";

$percent    = $infordb{'DNA'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'DNA'}->{'size'}/$total_repeat * 100);
print "Class II: DNA Transposon	$infordb{'DNA'}->{'num'}	$infordb{'DNA'}->{'size'}	$r_per	$percent\n";
my $tir_num = $infordb{'EnSpm/CACTA'}->{'num'} + $infordb{'hAT'}->{'num'} + $infordb{'MuDR/Mutator'}->{'num'} + $infordb{'TcMar'}->{'num'} + $infordb{'PIF-Harbinger'}->{'num'} + $otherDNA_num;
my $tir_size = $infordb{'EnSpm/CACTA'}->{'size'} + $infordb{'hAT'}->{'size'} + $infordb{'MuDR/Mutator'}->{'size'} + $infordb{'TcMar'}->{'size'} + $infordb{'PIF-Harbinger'}->{'size'} + $otherDNA_size;
my $percent = $tir_size / $g_size * 100;
$percent = sprintf("%.2f", $percent);
my $r_per = sprintf("%.2f", $tir_size / $total_repeat * 100);
print "         TIR     $tir_num        $tir_size       $r_per  $percent\n";
$percent    = $infordb{'EnSpm/CACTA'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'EnSpm/CACTA'}->{'size'}/$total_repeat * 100);
print "		EnSpm/CACTA	$infordb{'EnSpm/CACTA'}->{'num'}	$infordb{'EnSpm/CACTA'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'hAT'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'hAT'}->{'size'}/$total_repeat * 100);
print "		hAT	$infordb{'hAT'}->{'num'}	$infordb{'hAT'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'MuDR/Mutator'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'MuDR/Mutator'}->{'size'}/$total_repeat * 100);
print "		MuDR/Mutator	$infordb{'MuDR/Mutator'}->{'num'} $infordb{'MuDR/Mutator'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'TcMar'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'TcMar'}->{'size'}/$total_repeat * 100);
print "		Tc1/Mariner	$infordb{'TcMar'}->{'num'} $infordb{'TcMar'}->{'size'}	$r_per $percent\n";
$percent    = $infordb{'PIF-Harbinger'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'PIF-Harbinger'}->{'size'}/$total_repeat * 100);
print "		PIF/Harbinger	$infordb{'PIF-Harbinger'}->{'num'}	$infordb{'PIF-Harbinger'}->{'size'}	$r_per	$percent\n";
my $otherDNA_num  = $infordb{'DNA'}->{'num'} - $infordb{'EnSpm/CACTA'}->{'num'} - $infordb{'hAT'}->{'num'} - $infordb{'MuDR/Mutator'}->{'num'} - $infordb{'TcMar'}->{'num'}  - $infordb{'PIF-Harbinger'}->{'num'};
my $otherDNA_size = $infordb{'DNA'}->{'size'} - $infordb{'EnSpm/CACTA'}->{'size'} - $infordb{'hAT'}->{'size'} - $infordb{'MuDR/Mutator'}->{'size'} - $infordb{'TcMar'}->{'size'} - $infordb{'PIF-Harbinger'}->{'size'};
$percent    = $otherDNA_size/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$otherDNA_size/$total_repeat * 100);
print "		Other	$otherDNA_num $otherDNA_size	$r_per	$percent\n";
$percent    = $infordb{'Helitron'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'Helitron'}->{'size'}/$total_repeat * 100);
print "	Helitron	$infordb{'Helitron'}->{'num'}	$infordb{'Helitron'}->{'size'}	$r_per	$percent\n";
$percent    = $infordb{'tRep'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'tRep'}->{'size'}/$total_repeat * 100);
print "Tandem Repeats	$infordb{'tRep'}->{'num'}	$infordb{'tRep'}->{'size'}	$r_per $percent\n";
$percent    = $infordb{'unKnw'}->{'size'}/$g_size * 100;
$percent    = sprintf("%.2f",$percent);
$r_per      = sprintf("%.2f",$infordb{'unKnw'}->{'size'}/$total_repeat * 100);
print "Unkown	$infordb{'unKnw'}->{'num'}	$infordb{'unKnw'}->{'size'}	$r_per $percent\n"






