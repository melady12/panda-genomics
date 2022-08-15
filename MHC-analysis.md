```bash
seqkit grep -n -p Chr7 /mnt/analysis/genomicSource/v0-wild-panda/Ame_Qingling.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query mart_human.fasta -max_hsps 1 -outfmt 6 \
  | grep "HLA-" \
  | sort -k 9n \
  | wc -l
```

```bash
seqkit grep -n -p Chr7 /mnt/analysis/genomicSource/v0-wild-panda/Ame_Qingling.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query mart_mouse.fasta -max_hsps 1 -outfmt 6 \
  | grep "H2-" \
  | sort -k 9n > mouse_Qinling.blastn
 vi mouse_Qinling.blastn # counts
```

```bash
seqkit grep -n -p chr7 /mnt/analysis/genomicSource/v0-wild-panda/Ame_Sichuan.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query mart_mouse.fasta -max_hsps 1 -outfmt 6 \
  | grep "H2-" \
  | sort -k 9n > mouse_Sicuan.blastn
vi mouse_Sicuan.blastn  # counts
```

```bash
seqkit grep -n -p chr7 /mnt/analysis/genomicSource/v0-wild-panda/Ame_Sichuan.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query panda-MHC.fasta  -max_hsps 1 -outfmt 6 \
  |  sort -k 9n > panda_Sicuan.blastn
```

```bash
seqkit grep -p 5 /mnt/analysis/genomicSource/v2/Ailuropoda_melanoleuca.ASM200744v2.dna_sm.primary_assembly.all.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query mart_mouse.fasta -max_hsps 1 -outfmt 6 \
  | grep "H2-" \
  | sort -k 9n > mouse_v2.blastn
vi mouse_v2.blastn
```
```bash
seqkit grep -p 5 /mnt/analysis/genomicSource/v2/Ailuropoda_melanoleuca.ASM200744v2.dna_sm.primary_assembly.all.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query mart_human.fasta -max_hsps 1 -outfmt 6 \
  | grep "HLA-" \
  | sort -k 9n > human_v2.blastn
vi human_v2.blastn 
```
```bash
seqkit grep -p 5 /mnt/analysis/genomicSource/v2/Ailuropoda_melanoleuca.ASM200744v2.dna_sm.primary_assembly.all.fa \
  | blastn -subject - -subject_loc 20000000-35000000 -query panda-MHC.fasta -max_hsps 1 -outfmt 6  \
  | sort -k 9n > panda_v2.blastn
vi panda_v2.blastn 
```