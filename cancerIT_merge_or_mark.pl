# nice script but at the end of the day I don't know what options were used for my samples
# so this is not the most helpful thing I could have 

sub merge_or_mark_lanes {
  my ($options, @bams) = @_;
  my $tmp = $options->{'tmp'};

  my $marked = File::Spec->catdir($options->{'outdir'}, $options->{'sample'});
  
  if($options->{'cram'}) { $marked .= '.cram'; }
  else { $marked .= '.bam'; }

  return $marked if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my @commands = ('set -o pipefail');

  my $helper_threads = $options->{'threads'};

  my $input_str = join q{ }, sort @bams;

  my $strmd_tmp = File::Spec->catfile($tmp, 'strmdup');
  my $brc_tmp = File::Spec->catfile($tmp, 'brcTmp');

  my %tools;
  for my $tool(qw(bam_stats samtools md5sum)) {
    $tools{$tool} = _which($tool) || die "Unable to find '$tool' in path";
  }
  # md5sum is a programme that identifies 128-bit MD5 hashes 

  my $out_fmt = 'bam';
  my $idx_type = 'bai';
  my $idx_csi_flag = q{};
  if($options->{'cram'}) {
    $idx_type = 'crai';
    $out_fmt = 'cram';
    $out_fmt .= ',seqs_per_slice='.$options->{'seqslice'};
  }
  elsif(exists $options->{'csi'}) {
    # only valid for bam
    $idx_type = 'csi';
    $idx_csi_flag = '-c';
  }

  if(defined $options->{'nomarkdup'} && $options->{'nomarkdup'} == 1) {
      my $idx = q{};
      unless($options->{'noindex'}) {
        $idx = sprintf q{%s index -@ %d %s - %s.%s},
                        $tools{samtools}, $helper_threads, $idx_csi_flag, $marked, $idx_type;
      }

      my $namesrt = q{};
      $namesrt = q{-n} if($options->{'qnamesort'});

      my $merge    = sprintf q{%s merge %s -u -@ %d - %s},
                              $tools{samtools}, $namesrt, $helper_threads, $input_str;
      my $compress = sprintf q{%s view -T %s --output-fmt %s -@ %d -},
                              $tools{samtools}, $options->{reference}, $out_fmt, $helper_threads;
      my $md5      = sprintf q{%s -b > %s.md5},
                            $tools{md5sum}, $marked;
      my $stats    = sprintf q{%s -o %s.bas -@ %d},
                              $tools{bam_stats}, $marked, $helper_threads;
      push @commands, qq{$merge | pee "$stats" "$compress | pee '$idx' '$md5' 'cat > $marked'"};
  }
  else {
    my $merge;
    my $markdup;
    if(exists $options->{legacy}) {
      my $mmflagmod = _which('mmFlagModifier') || die "Unable to find 'mmFlagModifier' in path";
      my $bammarkdups = _which('bammarkduplicates2') || die "Unable to find 'bammarkduplicates2' in path";
      my $bammerge = _which('bammerge') || die "Unable to find 'bammerge' in path";
      my $bm_tmp = File::Spec->catfile($tmp, 'bmTmp');

      $input_str =~ s/ / I=/g;
      $merge = sprintf '%s SO=%s tmpfile=%s level=0 I=%s', $bammerge, 'coordinate', $bm_tmp, $input_str;

      my $mmQcRemove = sprintf '%s --remove -l 0 -@ %d', $mmflagmod, $helper_threads;
      my $bammarkdup = sprintf '%s tmpfile=%s M=%s.met level=0 markthreads=%d', $bammarkdups, $strmd_tmp, $marked, $helper_threads;
      my $mmQcReplace = sprintf '%s --replace -l 0 -@ %d', $mmflagmod, $helper_threads;

      $markdup = sprintf q{%s | %s | %s}, $mmQcRemove, $bammarkdup, $mmQcReplace;
    }
    else {
      $merge   = sprintf q{%s merge -u -@ %d - %s},
                         $tools{samtools}, $helper_threads, $input_str;
      $markdup = sprintf q{%s markdup --mode %s --output-fmt bam,level=0 -S --include-fails -T %s -@ %d -f %s.met - -},
                         $tools{samtools}, $options->{dupmode}, $strmd_tmp, $helper_threads, $marked;
    }
    my $compress = sprintf q{%s view -T %s --output-fmt %s -@ %d -},
                           $tools{samtools}, $options->{reference}, $out_fmt, $helper_threads;
    my $idx      = sprintf q{%s index -@ %d %s - %s.%s},
                           $tools{samtools}, $helper_threads, $idx_csi_flag, $marked, $idx_type;
    my $md5      = sprintf q{%s -b > %s.md5},
                           $tools{md5sum}, $marked;
    my $stats    = sprintf q{%s -o %s.bas -@ %d},
                           $tools{bam_stats}, $marked, $helper_threads;
    push @commands, qq{$merge | $markdup | pee "$compress | pee 'cat > $marked' '$idx' '$md5'" "$stats" };
  }

  if($options->{'cram'}) {
    push @commands, sprintf $CRAM_CHKSUM, $marked, $marked;
  }

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  return $marked;
}

# bammarkduplicates2 is likely refrring to this thing:
# https://manpages.debian.org/unstable/biobambam2/bammarkduplicates2.1.en.html 
# what this does: reads a coordinate sorted BAM file containing alignments computed by an aligner, 
# marks duplicate reads/alignments using the coordinates of the alignments and writes the marked alignments to a BAM file
# there aren't really that many options related to performance, so it shouldn't be that marking may be messed up by using some wrong settings 