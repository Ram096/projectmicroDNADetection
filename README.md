# MicroDNA Detection 

This project presents an algorithm for the detection of microDNAs using alignment data in Sequence Alignment/Map (SAM) format. By analyzing characteristics patterns in the alignment such as split, reads, soft-clipping and discordant mapping. The algorithm will identify genomic regions likely to produce circular DNA structures. 

## References 
Reference for the ([bam file])(https://drive.google.com/drive/folders/1WYPwCzSv__28iQlHwvZUZpNlmjTM7Fea).


## Algorithm

###### Read the raw bam file
$ python bamread.py
###### detect microDNAs
$ python detectmicro.py
###### plots
$ python plotmicroDNA.py



