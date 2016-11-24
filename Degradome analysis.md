Degradome softwares 

Tests and analysis

1. Introduction

[http://mirnablog.com/](http://mirnablog.com/)

[http://mirnablog.com/microrna-target-prediction-tools/](http://mirnablog.com/microrna-target-prediction-tools/)

2. Softwares

<table>
  <tr>
    <td>Soft</td>
    <td>Release date</td>
    <td>Installation</td>
    <td>Tested and working?</td>
    <td>pros / cons</td>
  </tr>
  <tr>
    <td>sPARTA</td>
    <td>2014</td>
    <td>Babel
/home</td>
    <td>yes</td>
    <td>Requires  Python3 (requires new environment)

Easy to use, fast (threading allowed)
The code is horrible...

Had to modify the code to accept all type of fasta…
Magic numbers everywhere in the script…</td>
  </tr>
  <tr>
    <td>Cleaveland</td>
    <td>2008</td>
    <td>Babel
/home</td>
    <td>yes</td>
    <td>(made by the same lab that developed ShortStack)
Livia doesn’t recommend this software</td>
  </tr>
  <tr>
    <td>PAREsnip
(srna-workbench</td>
    <td>2012</td>
    <td>MacOs
Babel
/home</td>
    <td>Works  on macOS
(with lot of rams)

Doesn’t work on babel...</td>
    <td>Requires last version of Java 1.8 or crash
Requires the Oracle java version, crashes sometime with the openjdk one…

I’ve installed all version of java on babel. From open jdk to oracle java version: type alternatives --config java and choose the version.
Even on the command line it requires a GUI to run...
Needs lots of memory (>16G)</td>
  </tr>
  <tr>
    <td>StarScan</td>
    <td></td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
</table>


3. Test 

1. **sPARTA** (**s**mall RNA-**PAR**E **T**arget **A**nalyzer)

**[https://github.com/atulkakrana/sPARTA.github/tree/master/spart**a](https://github.com/atulkakrana/sPARTA.github/tree/master/sparta)

Ce "soft" est juste un **script** à faire tourner sous **python3**. Du coup sur babel j’ai crée un environnement python3 (que j’ai appelé “**snakes**”) avec les packages nécessaire pour faire tourner le soft, ce qui me permet de garder toute mon install python2 et de  switch facilement en python3 quand il le faut.

Installation du nouvel environnement (snakes) sur babel:

**$ ****conda**** ****create**** --****name**** ****snakes python****=****3**** **

On switch:

**$ source activate snakes**

![image alt text](image_0.png)

Installation des packages de sPARTA:

$ conda install scipy

$ conda install -c bioconda pyfasta

$ conda install -c r rpy2

$ pip install rpy2 --upgrade

**USAGE**

python3 sPARTA.py -featureFile <featureFile.fa> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -libs <Lib_A.txt Lib_B.txt> -tarPred -tarScore --tag2FASTA --map2DD --validate

**GOOD TO KNOW**

* The filenames should **NOT** **contain** numbers

* The fasta files should be in **.fa**

* Every file must be in the same folder as the **sPARTA.py **script

* Every degradome libraries must be on** TAG COUNT** format (use tally)

* Need a folder per file. Can’t specify the output….

**# Tally command**

tally -i degradome.fastq -o degradome.unique.txt --nozip -format '%R%t%X%n'

* The miRNA.fa file header should not contain** "," **

This won’t work:

>sly-miR156d-3p,stu-miR156h-3p,stu-miR156i-3p,stu-miR156j-3p,stu-miR156k-3p

GCTCACTGCTCTATCTGTCACC

This will:

>sly-miR156d-3p

GCTCACTGCTCTATCTGTCACC

2. **Cleaveland**

[http://www.ncbi.nlm.nih.gov/pubmed/19017659](http://www.ncbi.nlm.nih.gov/pubmed/19017659)

A bit slow because it doesn’t have threading option but works. Generate a lot of T-plot, and a complete result file. 

Degradome alignment is made with Gstar and there is also the "mapping file"

$ CleaveLand4.pl -e degNbWT1.fasta -u nbe_nta_sly_stu_miRNA.fa -n Niben101_annotation.transcripts.fasta -o tutorial_plots > tutorial_results.txt

**GOOD TO KNOW**

Degradome files must be in **fasta format**

Gstar file:

T-plot:

![image alt text](image_1.jpg)

3. **PAREsnip**

**[http://srna-workbench.cmp.uea.ac.uk/doc/PAREsnip_UserGuide.pd**f](http://srna-workbench.cmp.uea.ac.uk/doc/PAREsnip_UserGuide.pdf)

**USAGE**

$ java -jar -Xmx4G ~/Softs/srna-workbenchV4.1.dAlpha/Workbench.jar -tool paresnip -srna_file sRNAseq/2_CleanReads/DrakarM.clean.fasta -deg_file Degradome/degDM.fa -tran_file ../ReferenceGenome/Brassica_napus/Index/Brassica_napus_v4.1.chromosomes.fa -genome_file ../ReferenceGenome/Brassica_napus/Brassica_napus_v4.1.chromosomes.fa -out_file test

The software should be used through the graphical user interface. The command line doesn’t work completely. It doesn’t produce the T-plots. 

**Command line with a config file: param.txt**

/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar -Xmx30g ~/Softs/srna-workbenchV4.2.1Alpha.D/Workbench.jar -tool paresnip -srna_file 18-30_sRNA/NbTMV_rep1_rep2.clean.fa -deg_file 18-30_degradome/NbTMV_rep1_rep2.clean.fasta -tran_file Niben101_annotation.transcripts.fasta -out_file TMV.paresnip.test.csv -params param.txt

