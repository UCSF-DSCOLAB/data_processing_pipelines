#!/usr/bin/perl -w

# code written by Ravi Patel

use strict;
use POSIX qw(pow log10 sqrt);
use FindBin;
use lib $FindBin::Bin;
use wGetOptions qw(wpod2usage wGetOptions);
use Time::HiRes;

my $vcf1 = "";
my $vcf2 = "";
my $field1 = "PL";
my $field2 = "PL";
my $gterror = 0.01;
my $region = "";

wGetOptions(
    "-Compare Identifies of Individuals from Inferred VCFs",
    "--Required Options",
    "vcf1=s" => [\$vcf1, "First Input VCF file. The FORMAT field must have PL or GT field. One of the VCF must have AF in the INFO field"],
    "vcf2=s" => [\$vcf2, "Second Input VCF file. The FORMAT field must have PL or GT field. One of the VCF must have AF in the INFO field"],
    "--Additional Options",
    "region=s" => [\$region, "Region to focus on"],
    "field1=s" => [\$field1, "FORMAT field entry to use as genotypes (GT) or genotype likelihoods (VCF)"],
    "field2=s" => [\$field2, "FORMAT field entry to use as genotypes (GT) or genotype likelihoods (VCF)"],
    "gt-error=f" => [\$gterror, "Assumed per-allele error rates When using GT field"]
    ) || wpod2usage(2);

unless ( ( $vcf1 ) && ( $vcf2 ) ) {
    print STDERR "ERROR: Missing required option\n";
    wpod2usage(2);
}

my %h = ();
my ($n1,$n2);
my $nv = 0;
my @id1s = ();
my @id2s = ();
my $errPL = int(0-log($gterror)/log(10));

die "ERROR: Cannot open $vcf1\n" unless ( -s $vcf1 );
die "ERROR: Cannot open $vcf2\n" unless ( -s $vcf2 );

if ( ( $vcf1 =~ /.vcf.gz$/ ) || ( $vcf1 =~ /.bcf$/ ) ) {
    open(V1,"bcftools view $vcf1 $region | grep -v ^##|") || die "Cannot open file\n";
}
else {
    open(V1,"cat $vcf1 | grep -v ^##|") || die "Cannot open file\n";
}
my @H1 = split(/[\t\r\n ]+/,<V1>);
@id1s = @H1[9..$#H1];
while(<V1>) {
    my ($chrom,$pos,$id,$ref,$alt,$qual,$filt,$info,$format,@G) = split;
    print STDERR "Reading $. variants at $chrom:$pos from $vcf1..\n" if ( $. % 10000 == 0 );             
    my $af = $1 if ( $info =~ /AF=([^;]+)/ );
    unless ( defined($n1) ) { $n1 = $#G+1; }
    my $key = "$chrom:$pos:$ref:$alt";
    my @PLs = ();
    my @formats = split(/:/,$format);
    my $ifield = -1;
    for($ifield=0; $ifield < @formats; ++$ifield) {
	last if ( $formats[$ifield] eq $field1 );
    }
    die "Cannot find field $field1 from $format\n" if ( $ifield > $#formats );
    for(my $i=0; $i < $n1; ++$i) {
	my @fields = split(/:/,$G[$i]);
	if ( $ifield <= $#fields ) {
	    if ( $field1 eq "PL" ) {
		my @p = split(/,/,$fields[$ifield]);
		push(@PLs,\@p);
	    }
	    elsif ( $field1 eq "GT" ) {
		my ($g1,$sep,$g2) = split(//,$fields[$ifield]);
		if ( $g1 eq "." ) { push(@PLs,[0,0,0]); }
		elsif ( $g1+$g2 == 0 ) { push(@PLs,[0,$errPL,2*$errPL]); }
		elsif ( $g1+$g2 == 1 ) { push(@PLs,[$errPL,0,$errPL]); }
		elsif ( $g1+$g2 == 2 ) { push(@PLs,[2*$errPL,$errPL,0]); }
		else { die "Cannot handle multi-allelic genotypes\n"; }
	    }
	}
	else {
	    push(@PLs,[0,0,0]);	    
	}
    }
    unless ( defined($h{$key}) ) {
	$h{$key} = [ undef, undef, $af ];
    }
    $h{$key}->[0] = \@PLs;
}
close V1;


if ( ( $vcf2 =~ /.vcf.gz$/ ) || ( $vcf2 =~ /.bcf$/ ) ) {
    open(V2,"bcftools view $vcf2 $region | grep -v ^##|") || die "Cannot open file\n";
}
else {
    open(V2,"cat $vcf2 | grep -v ^##|") || die "Cannot open file\n";
}
my @H2 = split(/[\t\r\n ]+/,<V2>);
@id2s = @H2[9..$#H2];
while(<V2>) {
    my ($chrom,$pos,$id,$ref,$alt,$qual,$filt,$info,$format,@G) = split;
    print STDERR "Reading $. variants at $chrom:$pos from $vcf2..\n" if ( $. % 10000 == 0 );     
    my $af = $1 if ( $info =~ /AF=([^;]+)/ );
    unless ( defined($n2) ) { $n2 = $#G+1; }
    my $key = "$chrom:$pos:$ref:$alt";
    my @PLs = ();
    my @formats = split(/:/,$format);
    my $ifield = -1;
    for($ifield=0; $ifield < @formats; ++$ifield) {
	last if ( $formats[$ifield] eq $field2 );
    }
    die "Cannot find field $field2 from $format\n" if ( $ifield > $#formats );
    for(my $i=0; $i < $n2; ++$i) {
	my @fields = split(/:/,$G[$i]);
	if ( $ifield <= $#fields ) {
	    if ( $field2 eq "PL" ) {
		my @p = split(/,/,$fields[$ifield]);
		push(@PLs,\@p);
	    }
	    elsif ( $field2 eq "GT" ) {
		my ($g1,$sep,$g2) = split(//,$fields[$ifield]);
		if ( $g1 eq "." ) { push(@PLs,[0,0,0]); }
		elsif ( $g1+$g2 == 0 ) { push(@PLs,[0,$errPL,2*$errPL]); }
		elsif ( $g1+$g2 == 1 ) { push(@PLs,[$errPL,0,$errPL]); }
		elsif ( $g1+$g2 == 2 ) { push(@PLs,[2*$errPL,$errPL,0]); }
		else { die "Cannot handle multi-allelic genotypes\n"; }
	    }
	}
	else {
	    push(@PLs,[0,0,0]);	    
	}
    }
    
    if ( defined($h{$key}) ) {
	$h{$key}->[1] = \@PLs;
	if ( defined($af) ) {
	    if ( defined($h{$key}->[2]) ) {
		$h{$key}->[2] = ($h{$key}->[2] + $af)/2;
	    }
	    else {
		$h{$key}->[2] = $af;
	    }
	}
	elsif ( !defined($h{$key}->[2]) ) {
	    die "ERROR: One of the VCFs must have AF field\n";
	}
    }
}
close V2;

### We need to calculate BF
my @llk0s = (0) x ($n1 * $n2);
my @llk1s = (0) x ($n1 * $n2);
my $l10 = log(10);

keys %h;
while( my ($k,$v) = each %h ) {
    my ($rPL1s, $rPL2s, $af) = @{$v};
    next unless ( defined($rPL1s) && defined($rPL2s) );
    next if ( ( $af == 0 ) || ( $af == 1 ) );
    my @PL1s = @{$rPL1s};
    my @PL2s = @{$rPL2s};

    my @gps = ( (1-$af)*(1-$af), 2*$af*(1-$af), $af*$af );
    my @p1s = ( $gps[0], 0, 0, 0, $gps[1], 0, 0, 0, $gps[2] );
    my @p0s = ( 0 ) x 9;
    for(my $i=0; $i < 3; ++$i) {
	for(my $j=0; $j < 3; ++$j) {
	    $p0s[$i*3 + $j] = $gps[$i] * $gps[$j];
	}
    }
    for(my $i=0; $i < $n1; ++$i) {
	for(my $j=0; $j < 3; ++$j) {
	    $PL1s[$i]->[$j] = exp(0-($PL1s[$i]->[$j])/10.0*$l10);
	}
    }
    for(my $i=0; $i < $n2; ++$i) {
	for(my $j=0; $j < 3; ++$j) {
	    $PL2s[$i]->[$j] = exp(0-($PL2s[$i]->[$j])/10.0*$l10);
	}
    }
    my @lk0s = (0) x ($n1 * $n2);
    my @lk1s = (0) x ($n1 * $n2);
    for(my $i=0; $i < $n1; ++$i) {
	for(my $j=0; $j < $n2; ++$j) {
	    for(my $g1=0; $g1 < 3; ++$g1) {
		for(my $g2=0; $g2 < 3; ++$g2) {
		    $lk0s[$i*$n2 + $j] += ( $p0s[$g1*3 + $g2] * $PL1s[$i]->[$g1] * $PL2s[$j]->[$g2] );
		    $lk1s[$i*$n2 + $j] += ( $p1s[$g1*3 + $g2] * $PL1s[$i]->[$g1] * $PL2s[$j]->[$g2] );		    
		}
	    }
	}
    }

    for(my $i=0; $i < $n1; ++$i) {
	for(my $j=0; $j < $n2; ++$j) {
	    $llk0s[$i*$n2 + $j] += log( $lk0s[$i*$n2 + $j] );
	    $llk1s[$i*$n2 + $j] += log( $lk1s[$i*$n2 + $j] );	    
	}
    }

    ++$nv;

    print STDERR "Calculating pairwise likelihoods for $nv variants..\n" if ( $nv % 10000 == 0 );
}

print STDERR "Finished processing $nv variants across ($n1,$n2) samples..\n";

print join("\t","SM_IDs",@id2s)."\n";
for(my $i=0; $i < $n1; ++$i) {
    print $id1s[$i];
    for(my $j=0; $j < $n2; ++$j) {
	#printf("\t(%.5lf,%.5lf)", $llk1s[$i*$n2 + $j], $llk0s[$i*$n2 + $j]);
	printf("\t%.2lf", $llk1s[$i*$n2 + $j] - $llk0s[$i*$n2 + $j]);	
    }
    print "\n";
}