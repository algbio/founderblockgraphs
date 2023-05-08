# xGFA format
**WARNING: extension under development, this is not a final version**

Extension of [GFA v1.1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for graphs derived from the segmentation of aligned sequences. Its intended use is to maintain the coordinates of the multiple sequence alignment in the graph (more specifically the starting columns of each segment).  

Example:
```
M	4	14
X	1	3	7	10
B	1	2	3	3
S	1	AG
S	2	CGA
S	3	CA
S	4	CTA
S	5	CTC
S	6	CT
S	7	GATAC
S	8	GTT
S	9	GTTAC
L	1	+	2	+	0M
L	1	+	3	+	0M
L	2	+	4	+	0M
L	2	+	5	+	0M
L	3	+	4	+	0M
L	3	+	6	+	0M
L	4	+	7	+	0M
L	4	+	8	+	0M
L	5	+	9	+	0M
L	6	+	9	+	0M
P	file1	1+,2+,4+,7+	*
P	file2	1+,3+,4+,8+	*
P	file3	1+,2+,5+,9+	*
P	file4	1+,3+,6+,9+	*
```

## `M` MSA line
Describes the multiple sequence alignment from which the graph was generated. `M	4	14` means that the alignment represented `4` sequences (rows) aligned into `14` positions (columns).

## `X` Segmentation line
Descibes the "cuts" of the MSA segmentation resulting in the graph, or more specifically the starting (1-indexed) column of each block. `X	1	3	7	10` means that the graph is made of 4 blocks, whose starting MSA columns are `1`, `3`, `7`, and `10` respectively.

## `B` Block lines
Describes the size of each block of the graph. `B	1	2	3	3` means that node `1` belongs to the first block, nodes `2` and `3` to the second block, nodes `4`, `5`, `6` to the third, and nodes `7`, `8`, `9` to the last.
