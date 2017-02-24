# samcompare
A tool to compare the output of multiple read mappers

This tool compares the output of different mappers on the same dataset, to report which tool reports the most reads as mapped/high quality and the difference in alignment.

##Licence

All parts of this tool is licenced under GPLv3.  
A copy of this licence is included under LICENSE.

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
