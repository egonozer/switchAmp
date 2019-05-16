# SwitchAmp

## 1. DESCRIPTION:  
A tool for characterizing antigenic switching in Neisseria gonorrhoeae pilus gene sequence from high throughput long-read (i.e. PacBio) amplicon sequences.

## 2. REQUIREMENTS: 
 
* Perl 5.10 or above
* Additional Perl modules 
  * Text::Levenshtein::XS
  * Inline (this is often part of the default installation on most systems)
	  
## 3. 	INSTALLATION:
Easy install:  
You may need to use `sudo` before `cpan` commands, depending on your setup.

	cpan install Text::Levenshtein::XS  
	cpan install Inline
	git clone ____________
	

## 4. USAGE:
### 4.1 Required Inputs:
`-q`  
Amplicon sequence reads in fastq format.  

* Read files can be gzipped as long as you have `gzip` installed on your system and in your PATH (this is default for most systems).  
* If you want to combine 2 or more read files, simply separate them with commas. For example:  
`-q read_file1.fastq,read_file2.fastq`  

`-n`  
File containing positions and sequences of non-conserved regions in silent and parental sequences. File format is as follows:  

* Header lines for each variable region start with `=` and are followed by these data, separated by commas:
  * Variable region number [integer]
  * ID of the parental sequence [string]
  * Position 1 base upstream of the variable region start, relative to the parental sequence. So, for example, if bases 5 through 9 in the parental sequence are variable, this value should be '4'. [integer]
  * Position 1 base downstream of the variable region stop, relative to the parental sequence. So, for example, if bases 5 through 9 in the parental sequence are variable, this value should be '10'. [integer]
  * Length of the variable region in the parental sequence, in bp [integer]
* Subsequent lines should include the ID of the parental or silent sequence and the variable region base or bases. ID and sequence should be separated by a tab character. If a variable region is missing, the sequence should be left blank. If a variable region represents a deletion relative to the parental, the sequence should be represented by "-" character.
* Example file:

```
=1,1-81-S2,166,169,2
1c1     GG
6c2     GG
6c3     
7c1     AG
1-81-S2 AG
=2,1-81-S2,222,235,12
1c1     TCCCCCGCCGACAAA
6c2     TCCCCCCCTCCGAC
6c3     TCCGCCTCCGAC
7c1     TCCTCCGCCACCGAC
1-81-S2 TCCGCTTCAACA
```
* This file can be automatically generated from a fasta-formatted aligment of parental and silent sequences using the included script `fasta_alignment_conserved.pl`

### 4.2 Options:
`-s`  
Parental sequence fasta file. This file should contain the full sequence of the parental gene corresponding to the coordinates given to the variable region file given to `-n`. If there are upstream and/or downstream sequences flanking the parental gene sequence in each amplicon sequence, these should also be included in this file with the headers ">up" and ">dn", respectively.

Default: The parental sequence of N. gonorrhoeae strain 1-81-S2 with fixed upstream and downstream sequences will be used. These sequences will be displayed by running the program with the `-S` (capital S) option.

`-S`  
Print the default parental, upstream, and downstream seqeuences and quit.

`-t`  
Maximum number of differences (mismactches + insertions + deletions) in upstream or downstream sequences.  
Default: 2

`-r`  
Minimum number of reads per unique sequence.  
Default: 3

`-d`  
Do not ignore variants in conserved regions.  
Default: Differences in conserved regions between the reference sequence and reads will be ignored.

`-p`  
Maximum number of differences in variable regions relative to the parental sequence for a read to be counted as parental. This will only apply if none of the silent copies are perfect matches in one or more of the variable regions.  
Default: 1

`-x`  
Maximum number of differences per variable region. If a read contains one or more variable regions that exceed this number of differences (mismatches + insertions + deletions) from all of the parental and silent sequences, the read will be omitted.  
Default: 2

`-v`  
File with list of silent sequenes associated with a phenotype of interest. Reads with variable regions containing one of these silent seqeunces will be marked in the output.

Format:

```
variant_number1:id1,id2,id3
variant_number2:id1,id3
```
* 'variant number' should match the numbering in the file given to `-n`
* 'id' should match sequence IDs in the file given to `-n`

`-o`  
Output file prefix.  
Default: "output"

`-h`  
Hierarchical clustering linkage method. Choices are 'single' or 'complete'.  
Default: complete

## 5. OUTPUTS:
* sequence_key.txt  
File containing all unique read sequences after removal of upstream and downstream sequences.

* results.only_variable.tsv  
File showing characterization of variable regions in each unique sequence.
  * Column "id": sequence ID. Corresponds to sequence IDs in "sequence_key.txt" file
  * Column "count": Number of reads represented by sequence
  * Column "has_var": Only present in file if option `-v` used. Shows sequences matching characteristics of interest.
  * Columns "varX:Y-Z": X is the variable region number, Y and Z are the beginning and end coordinates of the region relative to the parental sequence. Values in these rows show best matches with silent sequences. If blank, sequence in this region matched parental sequence with differences equal to or less than the threshold given by option `-p`. Otherwise, results will be formatted like this: "<# of differences>:<sequenceID1>,<sequenceID2>". For example, if the sequence of a variable region perfectly matches the silent sequence 2c1 only, the result will be "0:2c1". If the sequence of a variable region best matches both sequences 2c1 and 6c1 with 2 differences from each, the result will be "2:2c1,6c1". 
  * Column "top_vars": Lists the sequence or sequences matching the most number of variable regions. For example, if the sequence of region var1 matches 1c1 and 1c3 and the sequence of region var5 matches 1c1 and 2c2, the top_var will be listed as 1c1 because the sequence of two regions match this silent copy. 
  * Column "top_var_count": The number of variable regions matching the sequence or sequences in "top_vars". For the example above, the number would be 2.

* variants.only_variable.varX_Y-Z.fasta (multiple files)  
File containing all sequences not determined to be parental in a variable region where X is the variable region number, Y and Z are the beginning and end coordinates of the region relative to the parental sequence. The first sequence in the file will be the parental sequence. All subsequent sequences will have IDs corresponding to the unique sequence as found in the "sequence_key.txt" file and the first column of the "results.only_variable.tsv" file. 

* stats.txt  
File containing program run parameters and results of filtering and sorting steps.

## LICENSE:

SwitchAmp  
Copyright (C) 2019 Egon A. Ozer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  See LICENSE.txt

## CONTACT:

Contact [Egon Ozer](e-ozer@northwestern.edu) with questions or comments.
   
