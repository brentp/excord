usage
=====

By default excord expects stuff from a single chromosome:
```
excord --fasta $fasta --prefix some/cool/project $bam $chrom | bgzip -c > $out
```
This will write some/cool/project.rp.bin with reference and pair coverage interleaved.

It's possible to read a bam from stdin and write only the discorants and splitters:
```
<bam stream> | excord --discordantdistance 500 --fasta $fasta /dev/stdin | bgzip -c > $out
```
It's not possible to write reference coverage with this mode.
