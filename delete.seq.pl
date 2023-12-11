#!/usr/bin/perl
open FASTA, $ARGV[0];
$filename = (split /\./, $ARGV[0])[0];
open OUT, ">$filename.complete.seq";
while (<FASTA>){
chomp;
if (s/^>//){
@infor = split /\s/, $_;
$seqname = $infor[0]."_".$infor[1]."_".$infor[2]."_".$infor[3]."_".$infor[4];
}else {
$len{$seqname}  += length $_;
$fasta{$seqname} .= $_;
}
}
foreach $gene (keys %len){
if ($len{$gene} > 810){
print "$len{$gene}\n";
print OUT ">$gene\n$fasta{$gene}\n";
}
}
