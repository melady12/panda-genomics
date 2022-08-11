#!/usr/bin/perl
#==============================================================================#
#                                                                              #
#      Calculate contig and scaffold N50 to N90 value from sequence data       #
#          Version 1.70  by Zhang Zhouhao:zhangzhouhao@novogene.cn             #
#                Formated Output in Lines    2012-07-31                        #
#                                                                              #
#==============================================================================#

use 5.008;
use strict;
use warnings;
use Getopt::Std;
use PerlIO::gzip;

# usage
my $hello=<<HELLO;
Calculate contig and scaffold N50 to N90 value from sequence data (fasta/fa.gz format)
calSeq V1.7 2012-07-31 zhangzhouhao\@novogene.cn
Usage:
        perl $0 <Options>
Example:
        perl $0 -i Fasta.fa > result.xls
        perl $0 -i Fasta.fa.gz -o result.xls
Use -h option to show detailed options
HELLO

my $usage=<<USAGE;
calSeq V1.5 2012-05-08 zhangzhouhao\@novogene.cn
Usage:
        perl $0 <Options>
Options:	
	-i	(Req)	<String>	Diretory and file name of the input genome sequence.
					# Data file should be fasta style (.fa/.fasta/.fa.gz).
	-o	(Opt)	<String>	Diretory and file name prefix of the output genome sequence.
					# Default value is 'STDOUT'.
	-l	(Opt)	<int>		Cutoff of minium length for scaffold(length withoutN).
					# Default value is 100.
                                        # Now only filter for scaffold length
                                        # Calculate contig in scaffold after filter
        -c      (Opt)   <NONE>          Cut head and tail Ns of scaffold (default OFF)
        -t      (Opt)   <int[,int...]>  Output N[x]0 list, default (5,6,7,8,9){N50~N90}
	-h	(Opt)	<NONE>		Show detailed help message.
Example:
        perl $0 -i Fasta.fa > result.xls
        perl $0 -i Fasta.fa.gz -o result.xls
USAGE

# get options

my %opt;
getopts ('i:o:l:t:ch', \%opt);
die $hello if (keys %opt == 0);
die $usage if $opt{'h'};

# check options

die "No input file name!\n$usage" unless (defined $opt{'i'});
defined $opt{'l'} or $opt{'l'} = 100;
# warn $opt{'l'};
die "Invalid value of minium length cutoff!\n$usage" unless ($opt{'l'}>=0 and int($opt{'l'})==$opt{'l'});
defined $opt{'t'} or $opt{'t'} = "5,6,7,8,9";
my @n_lst = split ",",$opt{'t'};
foreach (@n_lst) {
        next if (/^[123456789]$/);
        die "Invalid value of N[x]s list!\n$usage";
}

# open handles

open IN,"<:gzip(autopop)",$opt{'i'} or die $!;
if (defined $opt{'o'}) {
        open OUT,">",$opt{'o'} or die $!;
        *STDOUT = *OUT;
}

# calculate length database from sequence

my @con;
my @sca;
my @lines;
my $con100 = 0;
my $sca100 = 0;
my $con2k = 0;
my $sca2k = 0;
my $totcon = 0;
my $totsca = 0;
while (<IN>) {
        chomp;
        next unless $_;
        if (/^>/) {
                if (@lines) {
                        my $curl = join '',@lines;
                        if ($opt{'c'}) {
                                $curl =~ s/^[Nn]*(.*)[Nn]*$/$1/;
                        }
                        @lines = ();
                        my $scal = length $curl;
                        if ($scal >= 100) {
                                ++$sca100;
                                if ($scal >= 2000) {
                                        ++$sca2k;
                                }
                        }
                        if ($scal >= $opt{'l'}) {
                                push @sca,$scal;
                                $totsca += $scal;
                        } else {
                                next;
                        }
                        my @ctgs = split /[Nn]+/,$curl;
                        foreach (@ctgs) {
                                my $conl = length $_;
                                push @con,$conl;
                                $totcon += $conl;
                                if ($conl >= 100) {
                                        ++$con100;
                                        if ($conl >= 2000) {
                                                ++$con2k;
                                        }
                                }
                        }
                }
        } else {
                push @lines,$_;
        }
}
my $curl = join '',@lines;
@lines = ();
my $scal = length $curl;
if ($scal >= 100) {
        ++$sca100;
        if ($scal >= 2000) {
                ++$sca2k;
        }
}
if ($scal >= $opt{'l'}) {
        push @sca,$scal;
        $totsca += $scal;
        my @ctgs = split /[Nn]+/,$curl;
        foreach (@ctgs) {
                my $conl = length $_;
                push @con,$conl;
                $totcon += $conl;
                if ($conl >= 100) {
                        ++$con100;
                        if ($conl >= 2000) {
                                ++$con2k;
                        }
                }
        }
        @ctgs = ();
}
close IN;

# sort length arrays

# die "no con\n" unless (@con);
my @con_s=sort {$b<=>$a} @con;
@con = ();
die "no sca\n" unless (@sca);
my @sca_s=sort {$b<=>$a} @sca;
@sca = ();

# print head text

my ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
($sec,$min,$hour,$mday,$mon,$year) = (
    sprintf("%02d", $sec),
    sprintf("%02d", $min),
    sprintf("%02d", $hour),
    sprintf("%02d", $mday),
    sprintf("%02d", $mon + 1),
    $year + 1900
);

my $headtext=<<HEADTEXT;
# Coverage of assemble result from file $opt{'i'}
# Data production: $year-$mon-$mday $hour:$min:$sec
# Minium length cutoff is $opt{'l'}
HEADTEXT
print $headtext;
my $out_title = sprintf "%-16s%-16s%-16s%-16s%-16s%-16s%-16s%-16s",'#Title','Total_length','Total_number','Num>=100','Num>=2000','Average_length','Max_length','Min_length';
for (@n_lst) {
        my $curTitleLen = 'N'.$_.'0_length';
        my $curTitleNum = 'N'.$_.'0_number';
        $out_title .= sprintf "%-16s%-16s",$curTitleLen,$curTitleNum;
}
$out_title =~ s/\s+$//;
print $out_title,"\n";

# calculate and print result

## Contig
my $res_con=&cal_cov(\@con_s,$totcon,\@n_lst);
my $out_contig = sprintf "%-16s%-16s%-16s%-16s%-16s%-16s%-16s%-16s",'Contig',${$res_con}{'tot_l'},${$res_con}{'tot_n'},$con100,$con2k,${$res_con}{'aver'},${$res_con}{'max'},${$res_con}{'min'};
for (@n_lst) {
        $out_contig .= sprintf "%-16s%-16s",${$res_con}{'l'.$_},${$res_con}{'n'.$_};
}
$out_contig =~ s/\s+$//;
print $out_contig,"\n";

## Scaffold
my $res_sca=&cal_cov(\@sca_s,$totsca,\@n_lst);
my $out_scaff = sprintf "%-16s%-16s%-16s%-16s%-16s%-16s%-16s%-16s",'Scaffold',${$res_sca}{'tot_l'},${$res_sca}{'tot_n'},$sca100,$sca2k,${$res_sca}{'aver'},${$res_sca}{'max'},${$res_sca}{'min'};
for (@n_lst) {
        $out_scaff .= sprintf "%-16s%-16s",${$res_sca}{'l'.$_},${$res_sca}{'n'.$_};
}
$out_scaff =~ s/\s+$//;
print $out_scaff,"\n";

if (defined $opt{'o'}) {
        close OUT;
}
# end and exit

exit 0;

# subroutines

sub cal_cov ($$$) {
	my ($arr,$len,$lst) = @_;
	my %cov;
	$cov{'tot_l'}=$len;
	$cov{'tot_n'}=@{$arr};
	$cov{'aver'}=int($cov{'tot_l'}/$cov{'tot_n'});
	foreach (@{$lst}) {
                ($cov{'l'.$_},$cov{'n'.$_})=&cov_n($arr,$cov{'tot_l'},$_);
        }
        $cov{'max'}=${$arr}[0];
        $cov{'min'}=${$arr}[-1];
        return (\%cov);
	exit 0;
}

sub cov_n ($$$) {
    my ($l_arr,$tot,$cov_p)=@_;
    my $limit=$tot*$cov_p/10;
    my $length=0;
    my $number=0;
    my $cur_len=0;
    foreach (@{$l_arr}) {
        last if ($length>=$limit);
        ++$number;
        $length+=$_;
        $cur_len=$_;
    }
    return ($cur_len,$number);
    exit 0;
}

__END__

============================
V1.70 2012-07-31
added column for shortest contig/scaffold output
