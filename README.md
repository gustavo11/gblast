# gblast

```
v0.7

gblast -q <FASTA of query sequences>  -r <m8 blast file> 

[ -s <FASTA of subject sequences> ] [-o <output png file>] [-w <width>] 

[-e <e-value threshold> ] [-c <coord. tolerance>] [ -g <group file>]


-o output file. default output name = <input file name>_[<group info>|<e_value>]_<coord. tolerance>.svg

-w width of final figure. default=1000.

-s FASTA of subject sequences. If not provided the figure will not show the 
length of subject sequence associated with each HSP.

-e quick way to set an overall cutoff for e-value. The e-value thresholds in 
the "group file" have precedence over this parameter.

-c <coord. tolerance>   . Indicates that HSP with higher scores should have 
priority during rendering.HSPs of lower score will not be rendered if they are <coord. tolerance> nt 
apart from HSPs with higher score. distance = max( start - start', end - end' )

Obs.: '-c 0' is different than omiting this parameter. '-c 0' means that HSPs of lower score will 
not be rendered if they overlap (distance = 0) with HSP of higher scores
If the '-c' parameter is omited, HSP's covering the exact same region will be shown.

-g <group file> defines group of subject sequences based on regex. E-value threshold and glyph 
color can be assigned to each group.
The definition file is a tab delimited file using the following format:

<group name>    [\w\W]+                 <color ex.:red,green,blue>      <e-value threshold ex.:0.1, 0.000001>
<group name>    gi\|74584587\|809998    <color ex.:red,green,blue>      <e-value threshold ex.:0.1, 0.000001>

Obs.: Please escape metacharacters in the regex, such as the '|' in the previous example

Output description:
-------------------
Description of items in the red caption associated to each HSP:
qc = coordinates in the query sequence.
sc = coordinates in the subject sequence.
sl = subject length.
s  = score.
l  = length of alignment in the query sequence.
e  = e-value.
%i = percentage of identity of the aligned fragment.
```
