#! user/bin/perl -w
use strict;

open (IN, $ARGV[0])||die;
my %hash;
my %partial;
my $scaffold = '';
my $Flag = 0;
my $title = '';

while (my $line = <IN>){
chomp $line;
       if ($line =~ /^>/){
            $scaffold = $line;
			$hash{$scaffold} = '';
            $Flag = 1;
       }
       elsif ($Flag == 1 && $line !~ /^>/){
               $hash {$scaffold} = $hash {$scaffold}. $line;
       }
       elsif ($Flag == 1 && $line =~ /^>/){
               $Flag = 0;
       }
}
close IN;

my @keys = keys %hash;
my $seq = '';
my @partial;
my @pseudo;

foreach my $key (@keys){
my $part_flag = 0;
my ($a, $b, $c) = split /\=/, $key;
my ($gi, $small, $large, $length) = split /\t/, $b;
my ($up, $down) = split /\t/, $c;
my $start = 0;
my $end = 0;
my $flag = 1;
if ($up < $down){
   $seq = $hash {$key};
   $hash {$key} = $seq;
   $start = $small;
   $end = $large;
}
if ($up > $down){
   $seq = $hash {$key};
   $seq = reverse $seq;
   $seq =~ tr/ATCGatcg/TAGCtagc/;
   $hash {$key} = $seq;
   $start = $length - $large + 1;
   $end = $start + $large - $small;
}
   my $t_start = $start + 27;
   while ($t_start < $end - 30){
   $t_start = $t_start + 3;
   my $code = substr ($seq, ($t_start - 1), 3);
      if ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i){
           $flag = 0;
           last;
      }
   }
   if ($flag == 1){
   my $s_flag = 0;
   my $e_flag = 0;
   my $n_start = -1;
   my $n_end = 5000;
   my $ns_flag = 0;
   my $ne_flag = 0;
           if ($start <= 30){
                $s_flag = 1;
           }
           if ($start > 30){
                my $s30 = $start;
                while ($s30 > $start - 31){
                $s30 = $s30 - 1;
                my $s30_code = substr ($seq, $s30, 1);
                    if ($s30_code =~ /N/i){
                         $n_start = $s30 + 1;
                         $ns_flag = 1;
                         last;
                    }
                }
           }
           if ($end >= $length - 30){
                $e_flag = 1;
           }
           if ($end < $length -30){
                my $e30 = $end - 1;
                while ($e30 < $end + 29){
                $e30 = $e30 + 1;
                my $e30_code = substr ($seq, $e30, 1);
                    if ($e30_code =~ /N/i){
                         $n_end = $e30 + 1;
                         $ne_flag = 1;
                         last;
                    }
                }
           }
           if ($s_flag == 1 || $ns_flag == 1 || $e_flag == 1 || $ne_flag == 1){
           my $START = $start + 30;
           my $END = $end - 33;
                if ($s_flag == 1 || $ns_flag == 1){
                    while ($START > 3 && $START > $n_start + 3){
                            $START = $START - 3;
                    }
                $start = $START;
                }
                if ( $s_flag == 0 && $ns_flag == 0){
                my $flaga = 0;
                $START = $start + 3;
                    while ($START > 3){
                    $START = $START - 3;
                    my $code = substr ($seq, ($START - 1), 3);
                        if ($code =~ /^ATG$/i){
                             $start = $START;
                             $flaga = 1;
                        }
                        elsif ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i || $code =~ /N/ig){
                             last;
                        }
                    }
                    if ($flaga == 0){
                         $start = $START;
                         $flaga = 1;
                    }
                }
                if ($e_flag == 1 || $ne_flag == 1){
                     while ($END <= $length - 6 && $END <= $n_end - 6){
                        $END = $END + 3;
                     }
                $end = $END + 3;
                }
                if ($e_flag == 0 && $ne_flag == 0){
                     my $flagb = 0;
                     while ($END <= $length - 6){
                     $END = $END + 3;
                     my $code = substr ($seq, $END, 3);
                         if ($code =~ /^TAG$/i || $code =~ /^TGA$/i || $code =~ /^TAA$/i || $code =~ /N/ig){
                              $end = $END + 3;
                              $flagb = 1;
                              last;
                         }
                     }
                     if ($flagb == 0){
                          $end = $END + 3;
                          $flagb = 1;
                     }
                }
           my $title = ">gi=$gi\t$start\t$end\t$length=$up\t$down";
           push @partial, $title;
           $partial {$title} = substr ($seq, ($start - 1), ($end - $start + 1));
		$part_flag = 1;
           }
   }
	if ($part_flag == 0){
		push @pseudo, $key;
		#$hash{$key} = substr ($seq, ($start - 1), ($end - $start + 1));
	}
}

open (PART, ">$ARGV[0]_partial.fas")||die;
open (PSEUDO, ">$ARGV[0]_pseudo_500extended.fas")||die;
foreach my $p_key (@partial){
        print PART $p_key, "\n", $partial{$p_key}, "\n";
}
close PART;

foreach my $pseudo (@pseudo){
	print PSEUDO $pseudo, "\n", $hash{$pseudo}, "\n";
}
close PSEUDO;


