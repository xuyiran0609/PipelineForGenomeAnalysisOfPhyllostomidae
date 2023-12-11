#!/usr/bin/perl
opendir A_SPECIES, $ARGV[0];
foreach $file (readdir A_SPECIES){
if ($file =~ /\.fas_pseudo_500extended\.fas$/ or $file =~ /\.fas_partial\.fas$/ or  $file =~ /\.complete.seq$/){
$filename = (split /\./, $file)[0];
open SEQ, "$ARGV[0]/$file";
open OUT,">$ARGV[0]/$filename.intact.complete.seq";
while (<SEQ>){
chomp;
if (s/^>//){
$name = $_;
}else {
$fasta{$name} .= uc ($_);
}
}
foreach $gene (keys %fasta) {
$start = substr ($fasta{$gene}, 0, 3);
$seq_reverse = reverse $fasta{$gene};
$end = substr ($seq_reverse, 0, 3);
if ($start =~ /ATG$/){
if ($end =~ /AAT$/ or $end =~ /GAT$/ or $end =~ /AGT$/){
print OUT ">$gene\n$fasta{$gene}\n";}
}
}
}
}
open OUTPEP, ">$ARGV[0]/$filename.intact.complete.pep.ok";
open OUTSH, ">$ARGV[0].transtation.sh";
print OUTSH "perl cds2aa.pl $ARGV[0]/$filename.intact.complete.seq > $ARGV[0]/$filename.intact.complete.pep\n"; 
system ("sh $ARGV[0].transtation.sh");
open PEP, "$ARGV[0]/$filename.intact.complete.pep";
while (<PEP>){
chomp;
if  (s/^>//){
$pepname = $ARGV[0]."_".(split /\s/, $_)[0];
}else {
$pep{$pepname} .= $_;
}
}
foreach $key (keys %pep){
if ($pep{$key} =~ /X/){
next;
}else{
print OUTPEP ">$key\n$pep{$key}\n";
}
}
