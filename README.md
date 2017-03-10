# samcompare
A tool to compare the output of multiple read mappers

This tool compares the output of different mappers on the same dataset, to report which tool reports the most reads as mapped/high quality and the difference in alignment.

##Licence

All parts of this tool is licenced under GPLv3.  
A copy of this licence is included under LICENSE.

##Dependencies

samcompare was developed using Python 3.5.2

##Help
samcompare is really simple to use:
The only thing needed is a list of name sorted sam files.

Example of usage:
```bash

python samcompare.py tool1.sam tool2.sam

```

Example on how to generate the sorted sam file:
```bash

tool1_mapping_command > tmp.sam
cat tmp.sam | elprep /dev/stdin /dev/stdout --sorting-order queryname --nr-of-threads 4 > tool1.sam
rm tmp.sam

python samcompare.py tool1.sam tool2.sam

```

###Parameters

positional arguments:
* multiple samfiles (at least 2)

optional arguments:
option | description
--- | ---
  -h, --help |      show this help message and exit
  -v, --verbose |
  -skip, --skip | Skip reads that are not found in all samples (if false, an error is occuring)
  -ccs, --ccs | These are PacBio CCS reads (fastq description should end with ccs (some mappers add extra numbers, these are removed to get the original query name))
  -rmpart RMPART | Remove this part of the file names (.sam is autoremoved)

