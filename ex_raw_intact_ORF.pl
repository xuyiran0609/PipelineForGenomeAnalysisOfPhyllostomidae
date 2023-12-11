#! user/bin/perl -w
use strict;

open (IN, $ARGV[0])||die;
my %hash;
my %genome;
my $scaffold = '';
my $flag = 0;
my $title = '';

while (my $line = <IN>){
chomp $line;
       if ($line =~ /^>/){
            $scaffold = $line;
			$hash{$scaffold} = '';
            $flag = 1;
       }
       elsif ($flag == 1 && $line !~ /^>/){
               $hash {$scaffold} = $hash {$scaffold}. $line;
       }
       elsif ($flag == 1 && $line =~ /^>/){
               $flag = 0;
       }
}
close IN;

my @keys = keys %hash;
my $seq = '';
my $state = 0;

open (OUT, ">$ARGV[0]_raw_intact.fas")||die;
open (OUU, ">$ARGV[0].no_intact.fas")||die;
foreach my $key (@keys){
$state ++;
#print "$state!\n$key\n";
        my ($a, $b, $c) = split /\=/, $key;
        my ($gi, $small, $large, $length) = split /\t/, $b;
        my ($up, $down) = split /\t/, $c;
        my $start = 0;
        my $end = 0;
        if ($up < $down){
             $seq = $hash {$key};
             $start = $small;
             $end = $large;
        }
        if ($up > $down){
              $seq = $hash {$key};
              $seq = reverse $seq;
              $seq =~ tr/ATCGatcg/TAGCtagc/;
              $start = $length - $large + 1;
              $end = $start + $large - $small;
        }
        my $STARTA = $start + 3;
        my $atg = 0;
        my $ttt = 0;
        my $flag1a = 0;
        my $flag1b = 0;
        while ($STARTA >= $start - 300 && $STARTA > 3){
        $STARTA = $STARTA -3;
                my $code = substr ($seq, ($STARTA - 1), 3);
                if ($code =~ /^ATG$/i){
                     $atg = $STARTA;
                     $flag1a = 1;
                }
                elsif ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i){
                        $flag1b = 1;
                        #$ttt = $STARTA + 3;
                        $ttt = $STARTA;
                        print "There found atg+\n";
                        last;
                }
        }
        my $STARTB = $start - 3;
        if ($flag1a == 0 && $flag1b == 1){
             #my $flag_t = 0;
             #my $startb = $STARTB;
             #while ($startb <= $start + 90){
             #$startb = $startb +3;
             #        my $code = substr ($seq, ($startb - 1), 3);
             #         if ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i){
             #             $flag_t = 1;
             #             $STARTB = $startb;
             #        }
             #}
             while ($STARTB <= $start + 90){
             $STARTB = $STARTB +3;
                     my $code = substr ($seq, ($STARTB - 1), 3);
                     if ($code =~ /^ATG$/i){
                          $atg = $STARTB;
                          $flag1a = 1;
                          last;
                     }
             }
             if ($flag1a == 0){
                  $atg = $ttt;
             }
        }
        elsif ($flag1a == 0 && $flag1b == 0){
                $atg = $STARTA;
        }
        my $ENDA = $end;
        my $stop = 0;
        my $flag2 = 0;
        while ($ENDA >= $end - 300){
        $ENDA = $ENDA -3;
        my $code = substr ($seq, $ENDA, 3);
            if ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i){
                 $stop = $ENDA + 3;
                 #print "There found stop coden-\n";
                 $flag2 = 1;
            }
        }
        my $ENDB = $end - 3;
        if ($flag2 == 1){
             #print "There found stop coden-\n";
        }
        elsif ($flag2 == 0){
             #while ($ENDB <= $end + 300 && $ENDB <= $length - 6){
		while ($ENDB <= $length - 6){
             $ENDB = $ENDB + 3;
                my $code = substr ($seq, $ENDB, 3);
                if ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i){
                     $stop = $ENDB + 3;
                     $flag2 = 1;
                     #print "There found stop coden-\n";
                     last;
                }
             }
		if ($flag2 == 0){
			$stop = $length;
			$flag2 = 1;
		}
        }
$title = "$a=$gi\t$atg\t$stop\t$length=$up\t$down";
$genome {$title} = substr ($seq, ($atg -1), ($stop - $atg +1));
}
my @g_keys = keys %genome;
foreach my $g_key (@g_keys){
my $g_flag = 1;
my $g_seq = $genome {$g_key};
my $g_len = length ($g_seq);
    if ($g_len < 750){
         $g_flag = 0;
    }
    elsif ($g_flag == 1){
            my $g_i = 0;
            while ($g_i < $g_len - 6){
            $g_i = $g_i + 3;
            my $g_code = substr ($g_seq, $g_i, 3);
                    if ($g_code =~ /^TAG$/i || $g_code =~ /^TGA$/i || $g_code =~ /^TAA$/i){
                         $g_flag = 0;
                    }
            }
    }
    if ($g_flag == 1){
         print OUT $g_key, "\n", $g_seq, "\n";
    }
    elsif ($g_flag == 0){
	my ($a, $b, $c) = split /\=/, $g_key;
	my ($gi, $d) = split /\t/, $b;
	foreach my $k (@keys){
		if ($k =~ /^>gi\=$gi.*\=$c$/){
			print OUU $k, "\n", $hash{$k}, "\n";
		}
	}
    }
}
close OUT;
close OUU;
