# This folder contains all the scripts used to run SnpEff

## Building the chetone database

The snpEff.config file was modified to include this line in the end:
```
Chetone_histrio.genome : Chetone histrio moth v1"
```

Then the database was build like this:

``` 
java -jar /mnt/scratch/projects/biol-specgen-2018/yacine/Tools/old_snpeff/snpEff.jar build -gff3 -v Chetone_histrio
```

