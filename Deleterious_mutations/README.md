# Building a SnpEff Database for a New Species
*(without CDS or protein FASTA files)*

This guide describes how to build a local SnpEff database for a species that does not have `cds.fa` or `protein.fa`.

---

## Requirements

- Java (required to run SnpEff)
- SnpEff
- samtools
- Reference genome FASTA
- GFF3 annotation file

---

## File naming conventions

Use the following filenames exactly:

- Reference genome: `sequences.fa`
- Annotation file: `genes.gff`

---

## Step-by-step instructions

### 1. Download the reference genome

Save the genome FASTA as:

```text
sequences.fa
```

---

### 2. Index the reference genome

```bash
samtools faidx sequences.fa
```

This step is mandatory.

---

### 3. Download the GFF3 annotation

Save the annotation file as:

```text
genes.gff
```

The GFF must correspond to the same assembly as `sequences.fa`.

---

### 4. Configure the SnpEff data directory

Edit `snpEff.config` and set:

```text
data.dir = /path/to/your/snpEff/data
```

---

### 5. Register the genome in `snpEff.config`

Add the following lines at the end of `snpEff.config`:

```text
# Genus species
Genus_species.genome : Genus_species v1
```

---

### 6. Create the genome directory

Inside `data.dir`, create the following directory:

```text
Genus_species/
```

Place these files inside it:

```text
sequences.fa
genes.gff
sequences.fa.fai
```

---

### 7. Build the SnpEff database

```bash
java -jar $SNPEFF build -gff3 -v Genus_species -noCheckCds -noCheckProtein
```

---

### 8. Verify the build

Ensure that the following file exists:

```text
Genus_species/snpEffectPredictor.bin
```

If this file exists and no fatal errors were reported, the database is ready to use.

---

## Notes

- `cds.fa` and `protein.fa` are **not required**.
- `-noCheckCds` and `-noCheckProtein` disable strict post-build validation.
- Java must be loaded before running SnpEff.
