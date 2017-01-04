# sam_variant_extractor
A small and specialised program to extract variants from sam files.

This is essentially the second version of my readSam program, but with a more
reasonable name, and uploaded to github for convenience. It is not my
intention to further develop this because:

1. It would be better to make use of htslib in order to read bam rather than
sam (primarily for space reasons)
2. The program is ridiculously specific (hence more of a script);
  1. it assumes and can only handle paired end sequence data
  2. the data must be sorted by the query name
  3. ? and probably has other constraints as well

The script makes use of memory mapping of the sam file to see if this allows
it to be speeded up using an openMP parallel for loop in the reading. However,
this doesn't seem to help at all; not sure why, but it seems that I may be
misunderstanding what mmap() does.
