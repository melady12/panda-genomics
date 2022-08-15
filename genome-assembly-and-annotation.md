### panda genome assembly pipeline tutorial
###  Step 1 定版基因组组装步骤、软件类型及软件参数
第一步: 组装contig序列
采用pacbio平台进行测序，测得三代数据量为282.45G，覆盖深度为113.89×（按照survey预估的基因组大小2.48G计算）。使用全部三代数据，利用falcon软件进行初步组装得到基因组的contigs序列。

```
#Title	Total_length	Total_number	Num>=2000	N50_length	N50_number	N90_length	N90_number
Contig	2,484,294,429	1,761	1,761	28,517,922	27	4,132,123	112
```
软件名：faclon
参数：
seed_coverage = 60
length_cutoff_pr = 7000
genome_size = 2500000000
pa_HPCdaligner_option =  -v -B128 -t16  -e.75  -k18 -h480  -l3000 -w8 -T8 -s1000
ovlp_HPCdaligner_option = -v -B128  -t16  -e.96 -k16  -h300  -l1000 -w8 -T16 -s1000

### Step 2 对基因组进行三代纠错
使用全部的三代数据，利用arrowr算法对基因组的contigs序列进行三代纠错，得到纠错后的基因组序列。

```
#Title	Total_length	Total_number	Num>=2000	N50_length	N50_number	N90_length	N90_number
Contig	2,488,096,155	1761	1761	28,564,904	27	4,137,256	112
```
算法：arrow
版本号：smrtlink_5.0.1
参数：默认参数

### Step 3 对基因组进行二代纠错
通过Illumina Hiseq测序平台进行双末端测序，获得的总测序量为450.86G，覆盖深度为181.80×。利用pilon软件，对基因组进行二代纠错。
```
#Title	Total_length	Total_number	Num>=2000	N50_length	N50_number	N90_length	N90_number
Contig	2,480,455,648	1761	1761	28,556,066	27	4,259,604	111
```
软件：pilon
版本号：pilon-1.22
参数：-Xmx300G --diploid --threads 20

### Step 4 连bionano
使用bionano策略平台获得总测序量313.15	G, 覆盖深度为126.27X
```
Total_length	Total_number	Max_length	N50_length	N50_number	N90_length	N90_number
Contig	2480475598	4610	101319274	28556066	27	4259604	111
Scaffold	2504077579	1585	193165166	58250622	16	11008006	50
```
软件名：BIONANO SOLVE
版本：Solve3.1
参数：
Denovo：-T 24 -j 5 -N 100 -i 3
hybrid-scaffold：-f -B 1 -N 1

### Step 5 连10X
使用10X Genomics策略平台获得总测序量293.19G, 覆盖深度为118.22X
```
Total_length	Total_number	Max_length	N50_length	N50_number	N90_length	N90_number
Contig	2480475598	4610	101319274	28556066	27	4259604	111
Scaffold	2504077579	1585	193165166	58250622	16	11008006	50
```
软件：fragScaff
版本：Version 140324.1
### Step 6 挂载染色体
根据测序得到的276.01  Hi-C数据，覆盖深度为111.29X使用Lachesis软件将组装得到的contigs/scaffolds序列提升到染色体水平，得到染色体基因组。
```
Total_length	Total_number	Max_length	N50_length	N50_number	N90_length	N90_number
Contig	2480475598	4654	101319274	28556066	27	4098111	111
Scaffold	2504105879	1382	202906273	134169173	8	43998035	18
```
软件：Lachesis
版本号：version-201701
RE_SITE_SEQ = GATC （酶）
CLUSTER_N = 21   （染色体数）
CLUSTER_MIN_RE_SITES = 200  （酶切位点数）
第七步：BAC评估
软件：LASTZ
参数：--seed=match12 --identity=85
