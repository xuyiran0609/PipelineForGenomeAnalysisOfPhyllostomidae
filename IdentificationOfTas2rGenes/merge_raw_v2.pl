#! user/bin/perl -w
use strict;

open (IN, $ARGV[0])||die;
my %hash;
my $state = 0;
my @gi;

while (my $line = <IN>){
chomp $line;
$state ++;
        $hash {$state} = $line;
         my ($a, $b, $c) = split /\t/, $line;
         #my ($d, $e, $f) = split /\|/, $b;
         push @gi, $b;
}
close IN;

my %count;
@gi = grep { ++$count{ $_ } < 2; } @gi;
my $gi_num = @gi;
my @keys = keys %hash;

my @temp;
foreach my $gi (@gi){
my @a_gi;
my @b_gi;
print "Here we are dealing with the gi $gi!\n";
        foreach my $key (@keys){
                my $value = $hash{$key};
                my ($A1, $B1, $C1, $D1, $E1, $F1, $G1, $H1, $START1, $END1, $Evalue1, $L1) = split /\t/, $value;
                #my ($d, $e, $f) = split /\|/, $B1;
                my $e = $B1;
                if ($e =~ /^$gi$/ && $START1 < $END1){
                        push @a_gi, $key;
                }
                elsif ($e =~ /^$gi$/ && $START1 > $END1){
                        push @b_gi, $key;
                }
        }
        my $a_num = @a_gi;
        my $b_num = @b_gi;
        print "$a_num + $b_num -\n";
        while ($a_num > 0){
        my @f_gi;
        print"1!\n";
                foreach my $gi_key (@a_gi){
                print "$a_num!\n";
                my $Value = $hash{$gi_key};
                my ($A, $B, $C, $D, $E, $F, $G, $H, $START, $END, $Evalue, $L) = split /\t/, $Value;
                my $EVALUE = $Evalue;
                my $LENGTH = $H - $G + 1;
                @f_gi = grep (!/^$gi_key$/, @a_gi);
                @a_gi = @f_gi;
                my $small = $START;
                my $large = $END;
                my $up = $START - 1000;
                my $down = $END + 1000;
                           my @A_gi = @a_gi;
                           foreach my $temp_gi_key (@A_gi){
                           my $temp_value = $hash{$temp_gi_key};
                           my ($a, $b, $c, $d, $e, $f, $g, $h, $start, $end, $evalue, $l) = split /\t/, $temp_value;
                           my $length = $h - $g + 1;
                              if ($start > $up && $end < $down && $start < $end){
                              @f_gi = grep (!/^$temp_gi_key$/, @a_gi);
                              @a_gi = @f_gi;
                                      if ($evalue < $EVALUE){
                                               $small = $start;
                                               $large = $end;
                                               $EVALUE = $evalue;
                                               $LENGTH = $length;
                                         }
                                      elsif ($evalue == $EVALUE && $length > $LENGTH){
                                              $small = $start;
                                              $large = $end;
                                              $EVALUE = $evalue;
                                              $LENGTH = $length;
                                      }
                              }
                           }
                my $temp = "$gi\t$small\t$large\t$EVALUE";
                push @temp, $temp;
                $a_num = @a_gi;
                print "$a_num!!\n";
                }
        $a_num = @a_gi;
        }
        while ($b_num > 0){
        my @f_gi;
        print "2!\n";
                 foreach my $gi_key (@b_gi){
                 print "$b_num!\n";
                 my $Value = $hash{$gi_key};
                 my ($A, $B, $C, $D, $E, $F, $G, $H, $START, $END, $Evalue, $L) = split /\t/, $Value;
                 my $LENGTH = $H - $G + 1;
                 my $EVALUE = $Evalue;
                 @f_gi = grep (!/^$gi_key/, @b_gi);
                 @b_gi = @f_gi;
                 my $small = $END;
                 my $large = $START;
                 my $up = $START + 1000;
                 my $down = $END - 1000;
                 my @B_gi = @b_gi;
                          foreach my $temp_gi_key (@B_gi){
                          my $temp_value = $hash{$temp_gi_key};
                          my ($a, $b, $c, $d, $e, $f, $g, $h, $start, $end, $evalue, $l) = split /\t/, $temp_value;
                          my $length = $h - $g + 1;
                               if ($start < $up && $end > $down && $start > $end){
                                       @f_gi = grep (!/^$temp_gi_key$/, @b_gi);
                                       @b_gi = @f_gi;
                                       if ($evalue < $EVALUE){
                                                $small = $end;
                                                $large = $start;
                                                $EVALUE = $evalue;
                                                $LENGTH = $length;
                                                }
                                       elsif ($evalue == $EVALUE && $length > $LENGTH){
                                               $small = $end;
                                               $large = $start;
                                               $EVALUE = $evalue;
                                               $LENGTH = $length;
                                               }
                               }
                          }
                 my $temp = "$gi\t$large\t$small\t$EVALUE";
                 push @temp, $temp;
                 $b_num = @b_gi;
                 print "$b_num!\n";
                 }
        $b_num = @b_gi;
        }
}

my %number;
@temp = grep { ++$number{ $_ } < 2; } @temp;
#print "@temp \n";

open (OUT, ">$ARGV[0]_merged")||die;
foreach my $out (@temp){
        print OUT $out, "\n";
}
close OUT;


