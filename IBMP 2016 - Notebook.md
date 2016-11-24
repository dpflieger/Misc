IBMP 2016 - Notebook

### Mardi 01 mars 2016

Configuration MacOs + installation des softs 

Paperasses CNRS (--> visite médicale à faire)

Biblio sur les smallRNAseq (Mendeley + Google doc)

Code photocopieur: 	287 394

Connexion mac:	ssh [dpflieger@insilico.u-strasbg.fr](mailto:dpflieger@insilico.u-strasbg.fr)

VPN: 			vpn.u-strasbg.fr

**Database plantes :**

TAIR10 		[https://www.arabidopsis.org](https://www.arabidopsis.org)

Phytozome		[https://phytozome.jgi.doe.gov/pz/portal.html](https://phytozome.jgi.doe.gov/pz/portal.html)

ENSEMBL plants	[http://plants.ensembl.org/index.html](http://plants.ensembl.org/index.html)

miRNA			[http://www.mirbase.org](http://www.mirbase.org)

**Liens intranet IBMP :	**

[https://ibmp-intranet.u-strasbg.fr/wiki/projects/intranet/Intranet.html](https://ibmp-intranet.u-strasbg.fr/wiki/projects/intranet/Intranet.html)

[https://ibmp-intranet.u-strasbg.fr/wiki/projects/ibmpwebapp/IBMP_WebApps.html](https://ibmp-intranet.u-strasbg.fr/wiki/projects/ibmpwebapp/IBMP_WebApps.html)

### Mercredi 02 mars 2016

Dicer-like protein, in plant there are 4 types going from DCL1 to DCL4. They can process dsRNA and produce small interfering RNA in different size.

![image alt text](image_0.png)

![image alt text](image_1.png)

Exemple de hairpin des pré-mir

![image alt text](image_2.png)

![image alt text](image_3.png)

Presentation about miRNAs and degradome sequencing:

[http://www.lcsciences.com/documents/degradome-sequencing-webinar.pdf](http://www.lcsciences.com/documents/degradome-sequencing-webinar.pdf)

### Jeudi 03 mars 2016

Prise de connaissance du project LabCom entre [Plant Advanced Technologies](http://www.plantadvanced.com/) et l’IBMP.

Intérêt pharmaceutique et cosmétique des terpènes et terpénoïdes. Détient le brevet PAT (Plantes à Traire)

Personnes impliquées:

* Danièle Werck-Reichhart

* Nicolas Navrot 

* Frédéric Bourgaud 

* Jean-François Ginglinger

* Claire Parage (Post-doc)		

Pour le projet → séquençage d’ARN de Scrophularia Nodosa avec et sans traitement hormonale (x3 réplicats) →  RNAseq paired-end, 125 bp en HiSeq-2500 séquencé chez Eurofins. 

Valérie m’a envoyé les code d’accès pour le site 1kp.

		 	

Mise en place de mes repos git pour l’IBMP.

Listage des outils de détection et d’analyse des small RNAs

* Info sur les [réseaux métabolomiques](https://en.wikipedia.org/wiki/Metabolic_network_modelling)

Listage des databases de plantes

Link about genome automatic annotation:	[http://www.ncbi.nlm.nih.gov/books/NBK20253](http://www.ncbi.nlm.nih.gov/books/NBK20253/)

### Vendredi 04 mars 2016

(Réunion avec Claire Parage pour le LabCom --> décalée)

Installation sur l’ordi à Valérie. Re-configuration et re-mise en place des softs, scripts, etc…

Ré-installation d’anaconda et tout le tralala pour python/perl.

Initiation à Gbrowse2 [https://bib.oxfordjournals.org/content/early/2013/02/01/bib.bbt001.full](https://bib.oxfordjournals.org/content/early/2013/02/01/bib.bbt001.full)

Récupération des fichiers gff3 sur le ftp de TAIR10. Il y a également un fichier de config pour Gbrowse2. → On pourra s’en servir comme exemple.

Pour le Gbrowser, il faut mettre en place une base mysql ou sqlite car la mémoire ne suffit pas pour load** TAIR10 **à la volée (ça rame à mort)

Check les browsers de l’UCSC pour voir ce qui est faisable.

### Lundi 07 mars 2016

Accès à la dropbox IBMP → afp:[//dpflieger@dropbox.u-strasbg.fr](mailto://dpflieger@dropbox.u-strasbg.fr)

Tutorial fait sur GBrowse2 + création de la base mysql sur la VM

Réunion avec Valérie sur les différents projets:

**I. Projet de Manfred**: Virus transcriptional impact in 2 differents plant species (mostly silencing)

Ils ont fait du RNAseq, small RNAseq, degradome. (Pas de réplicat...)

→ Regarder avec Stéfanie le 18 mars (projet en suspens pour l’instant)

→ Regarder aussi avec le doctorant Nicolas Pitzalis.

Faire une doc sur l’analyse du dégradome, les outils utilisés, tout ça…

**II. Projet de Claire Parage**: PAT annotation de pathways, optimisation de production de certaines molécules sur Scrophularia nodosa.

**III. Projet Annotation des TRFs** → download des banques (Valérie en a déja téléchargé trois sur son ordi)

**IV. Projet Cyrille small RNA **→ pour l’instant Timothée s’en est occupé

+ Projet annexe : configuration de GBrowse2

### Mardi 08 mars 2016

Avec Nicolas Pitzalis, nous avons fait un point de son projet sur l’étude des virus hôtes TuMV chez le colza (*Brassica napus*) et TMV chez la plante du tabac (*Nicotiana benthamiana)*

Ils ont fait du RNAseq, small RNAseq, ainsi que du dégradome avec un laboratoire à Madrid.

La boîte qui a séquencée les samples est [Ascidea](http://www.ascidea.com/). 

Le but de son PhD est de trouver des gènes candidats nécessaire au virus pour se propager et survivre ainsi que les gènes/pathways de la plante impliqués dans la défense contre ce virus. 

→ The virus-encoded suppressors of RNA silencing target

**La réference pour ****_Arabidopsis thaliana_**** est Columbia-0 (Col-0) **

Small RNA Truseq Adapters from Illumina documentation

RNA 5’ Adapter (RA5), part # 15013205

5’ GUUCAGAGUUCUACAGUCCGACGAUC

RNA 3’ Adapter (RA3), part # 15013207

5’ TGGAATTCTCGGGTGCCAAGG

Si jamais les outils pour l’analyse du dégradome et de small RNA sont décrits sur le site : [OmicTools](http://omictools.com/)

"The *Brassica* species include an important group of vegetable and oil crops and their genomes have complex evolutionary histories. A major focus for research has been **_Brassica napus_** (oilseed rape). This is an allopolyploid species formed by the hybridization of progenitor species **_Brassica rapa_** (which contributed the A genome) and **_Brassica oleracea_** (which contributed the C genome)"

Il faut que je dl le génome de référence de Brassica napus ainsi que celui de Nicotiana benthamiana.

### Mercredi 09 mars 2016

Pas d’installation de XtraFinder pour l’instant sur ElCaptain… dommage.

Du coup, tweak MacOS pour la commande "Copier le chemin" de Xtrafinder: [https://retinaboys.com/2015/07/19/comment-copier-chemin-fichier-dossier-osx/](https://retinaboys.com/2015/07/19/comment-copier-chemin-fichier-dossier-osx/)

Installation de [ZSH](https://github.com/robbyrussell/oh-my-zsh) et des plugins MacOS.

**Litterature:**

**Garcia-Ruiz H, Carbonell A, Hoyer JS, Fahlgren N, Gilbert KB, Takeda A, et al. (2015) Roles and Programming of Arabidopsis ARGONAUTE Proteins during Turnip Mosaic Virus Infection. PLoS Pathog 11(3): e1004755. doi:10.1371/journal. ppat.1004755**

"RNA silencing is a primary, adaptive defense system against viruses in plants. Viruses have evolved counter-defensive mechanisms that inhibit RNA silencing through the activity of silencing suppressor proteins. Understanding how antiviral silencing is controlled, and how suppressor proteins function, is essential for understanding how plants normally resist viruses, why some viruses are highly virulent in different hosts, and how sustainable antiviral resistance strategies can be deployed in agricultural settings"

**Genome:**

Pour le génome de référence de *Brassica napus*: assemblé et annoté par le genoscope, dispo ici :	[http://www.genoscope.cns.fr/brassicanapus/](http://www.genoscope.cns.fr/brassicanapus/)

Pour le génome de référence de *Nicotiana benthamiana *:

[http://sefapps02.qut.edu.au/benWeb/subpages/downloads.php](http://sefapps02.qut.edu.au/benWeb/subpages/downloads.php)

ftp://[ftp.solgenomics.net/genomes/Nicotiana_benthamiana](http://ftp.solgenomics.net/genomes/Nicotiana_benthamiana)

Mais à vérifier avec Stéfanie, pour l’instant je suis bloqué niveau data… Je n’ai accès à rien du tout.

Test de Abyss-Trans + Trinity en attendant.

### Jeudi 10 mars 2016

Accès au données de Stéfanie → checking des pipelines d’analyse

Regarder le sofware "patman" 

 

Pour installer la banque nr/nt/othergenomic etc... check ici → [http://www.ncbi.nlm.nih.gov/public/](http://www.ncbi.nlm.nih.gov/public/)

**Récap discussion avec Claire**:

**TerpFactory project** (un labcom entre PAT et l’IBMP)

Ils recherchent la voie de synthèse exact (gènes impliqués) chez *Scrophularia nodosa** *capable de synthétiser la molécule "[Harpagoside](http://www.scbt.com/datasheet-203073-harpagoside.html)" (un iridoid glycosides, anti-inflammatoire, antibactérien) prisé par les industries cosmétiques et pharmaceutiques. Un pathway est déjà connu et publié pour une molécule semblable chez *Catharantus roseus.*

Du coup, elle a séquencée 2 échantillons de RNA (avec 3 réplicats) sous 2 conditions (induced or non-induced by MeJa):

1. MeJa 10% in EtOH (methyl jasmonate)

2. EtOH 100%

(voir devis eurofins pour + d’info)

The gene in the harpagoside are downregulated, so normally we should observe a negative regulation of gene expression after treatment with MeJa.

WCBD:

1. Comparaison de la structure de la molécule finale avec d’autres molécules similaires dont les pathways sont connus.

2. Annotation & Recherche de similarité des gènes dans les banques. 

3. Expression & Co-expression des gènes.

### Vendredi 11 mars 2016

Installation de **blast+ **sur ma machine ainsi que **[sequenceserve**r](http://www.sequenceserver.com/).

EDIT: J’ai également installer **sequenceserver** sur l’ordi de Claire et celui de Nicolas. Ils pourront blaster sur les assemblages* de novo* du coup. Voir avec Francois pour avoir une VM avec **sequenceserver**, comme pour Gbrowse et JBrowse.

Downloaded:

* banque NR/NT  (la commande pour update:  update_blastdb.pl nr)

* pfam-A

* Uniref90

* Rfam

* uniprot_swissprot

* PhytozomeV11

* miRBase

Pour construire les banques avec **blast+**  → **makeblastdb**

Outils à test si jamais:

**Trapid:**

[http://bioinformatics.psb.ugent.be/webtools/trapid/](http://bioinformatics.psb.ugent.be/webtools/trapid/)

**Trinotate:**

[https://trinotate.github.io/](https://trinotate.github.io/) 

**Trinity guide on RNAseq : **

[http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/)

[https://github.com/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/wiki/Trinity-De-novo-Transcriptome-Assembly-Workshop](https://github.com/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/wiki/Trinity-De-novo-Transcriptome-Assembly-Workshop)

ftp://[ftp.broadinstitute.org/pub/users/bhaas/rnaseq_workshop/rnaseq_workshop_2014/Trinity_workshop_activities.pdf](http://ftp.broadinstitute.org/pub/users/bhaas/rnaseq_workshop/rnaseq_workshop_2014/Trinity_workshop_activities.pdf)

**RNAseq annotation :**

[https://www.biostars.org/p/108335/](https://www.biostars.org/p/108335/) 

Dommage pas de génome de référence pour *Scrophularia nodosa* (nb chromo? genome size?)

### Lundi 14 mars 2016

Check du package [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#abstract) pour l’analyse de co-expression des gènes

Installation Deseq2 + EdgeR + limma + CummRbund

Trinity fonctionne sur babel, installation de Trinotate

Install des banques également sur babel

Les banques ont été installé sur **/mnt/data/bioimage/databases. **

(par contre pas d’accès quand on lance on job en sbatch, faudra les mettre autre part)

Check le pipeline "Eugene" de l’INRA pour faire de l’annotation automatique après assemblage.

A l’INRA, ils utilisent également Artemis, mais ça a l’air plutôt archaîque.

Thèses sur l’étude de TuMV / TMV chez Thaliana et Nicotiana

[http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=4056&context=etd](http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=4056&context=etd)

[http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=16587&context=rtd](http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=16587&context=rtd)

→ Pas mal d’info intéressantes

### Mardi 15 mars 2016

Biblio

### Mercredi 16 mars 2016

Analyse des données d’Esther Lechner:

Map des 2 listes des gènes HB16 et HB6 sur la banque InterMine de phytozome ainsi que sur leur bioMart et récupération des informations.

Synthèse sur le projet avec Claire pour la réunion de demain

### Jeudi 17 mars 2016

**Réunion avec Claire et Valérie sur le projet TerpFactory**

→ discuter avec Eurofins pour l’assemblage de novo des transcripts

**Discussion avec Valérie sur son projet avec les ****[tRF**s](http://drive.google.com/open?id=1q-mBOry7jMuXFscra0_DmgTJFWiUdTf6318kQlRLDjo)

Download de la banque MPGR dans **/databases**

Installation des libs pour le support 32 bit des logiciels sur Babel

### Vendredi 18 mars 2016

**Réunion avec Manfred / Stéfanie / Nico / Khalid. **

**Projet assez complexe, pas vraiment de punchline, données réparties un peu partout,**

**et pas beaucoup de résultats positifs pour l’instant. **

**J’ai reçu aucun document à part celui de Livia et de Stéf…**

![image alt text](image_4.png)

### Lundi 21 mars 2016

J’ai appelé Nathalie Potier d’Eurofins pour discuter des assemblages de novo du projet de Claire. Après discussion avec la nana, j’ai écris à l’équipe de bioinfo en charge de l’analyse pour voir s’ils sont d’accord changer les modalités d’assemblage. 

Du coup,

1. MeJa (pool des 3 réplicats)

2. NT (pool des 3 réplicats)

3. MeJa_NT (pool des 6 samples)

Les assemblages vont être fait avec Velvet + Oases.

Check la publication sur l’analyse des 52 *Brassica napus *

[http://www.nature.com/articles/sdata201572](http://www.nature.com/articles/sdata201572)

Il y a les données sous forme de matrice, ainsi que les fichiers vcf correspondants pour chaque souche.

### Mardi 22 mars 2016

J’ai lancé le mapping des RNAseq de Manfred:

* avec 1 mismatch

* 2 mismatch

* 1 mistmatch nomultihits 

* Sans gff en input pour voir

* Avec les mêmes paramètres de Stéf pour comparer

* 3 mismatches

Ne sachant pas si les données étaient généré avec un kit Illumina "strand-specific" j’ai également fait plein de tests. 

EDIT: pas strand-specific

→ J’ai lancé les assemblages de novo des unmapped avec SOAPdenovo pour voir également ce qui peut sortir

EDIT: plein de mito / plein de chloro. Ca doit être normal pour les plantes → à voir !

Installation de brew sur le mac, puis install de samtools et tabix, et de la commande tree:

brew tap homebrew/science

brew install tabix

brew install samtools

brew install tree

**Mercredi 23 mars 2016 & Jeudi 24 mars 2016**

Checking des mapping des RNAseq

Mise en place du pipeline GATK pour le SNP calling

I made an EndNote file with a full report on mapping parameters tested and mapping stats

**Mardi**** 29 mars 2016**

Grosse discussion avec Nicolas P. sur *Brassica napus & Nicotiana Benthamiana*

**Mercredi 30 mars 2016**

For the TruSeq Small RNA Sample Prep Kit the adapter sequence is 5’ **TGGAATTCTCGGGTGCCAAGG**

**Illumina doc on adapters:**

**[http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pd**f](http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf)** **

**Attention**: adapteurs présents dans les fichiers clean de Stéfanie, ainsi que dans les fichiers clean de Livia…

Reclean des données de small + mapping sur TuMV 0 mismatch + mapping Brassica 0 mismatch + Stats mapping

Mapping des RNAseq sur le génome de ZS11, (génome de référence apparemment utilisé par l’équipe en Espagne...) par contre pas annoté du tout… erf....

**Jeudi 31 mars 2016**

**Vendredi 1 avril 2016**

Install et test de SOAP2aligner: impossible à faire tourner correctement →  every fu***ing time une segmentation fault...

**Lundi 4 avril 2016**

Outils pour analyser les miR 

[http://crdd.osdd.net/rna.php](http://crdd.osdd.net/rna.php)

**Mardi 4 avril 2016**

Analyse du RNAseq de *Brassica napus* avec EdgeR, ainsi qu’avec DEseq2.

Installation de PareSNIP sur mon ordi pour refaire l’analyse des targets et du dégradome.

Livia fini son post-doc bientôt et je pense que Manfred aimerait bien qu’on sache faire l’analyse en local ici.

TODO: Regarder et faire un point sur : Trinotate pipeline, RSEM, database

Pour Manfred:

Tout compte fait, il aimerait utiliser les transcriptomes **exactes** de Tanto et Drakar pour la suite de l’analyse. L’idéal aurait été de faire de l’assemblage *de novo* dès le début (je vais le faire en parallèle si j’ai le temps au cas où), mais bon…

Du coup, j’ai mappé les transcripts Tanto/Drakar sur le transcriptome de la référence Darmor-bzh, puis j’ai call les SNPs le pipeline GATK-RNSeq. J’ai ensuite inféré les SNPs identifiés dans le transcriptome de Darmor-bzh, pour avoir un Drakar-like transcriptome et Tanto-like transcriptome. Vu que l’on utilise du RNAseq (en plus sans réplicats), les SNPs called ne sont pas fiables du tout. J’ai d’ailleurs comparé mes positions polymorphiques avec ceux du papier des 52 brassica napus. Seulement ~10k SNPs sont en correspondance entre tous les cultivars…

Dans le papier, ils ont comparé leurs SNPs avec le 60k SNP Brassica Consortium Infinium genotyping:

TODO → essayer de récup les données 60k (si possible) et faire aussi une comparaison.

EDIT → impossible de trouver les données en ligne… ce n’est pas encore publié. Il faut apparemment faire une demande au consortium…

Après re-discussion avec Manfred, il va séquencer le full génome de Tanto et Drakar. (ce qui n’est pas plus mal si Nico veut une analyse des SNPs/Indel).

**Mercredi 5 avril 2016**

Réception des données de séquençage de Claire sur DD externe. Sauvegarde sur mon ordi ainsi que sur babel :

/mnt/data/bioimage/NGS-archives/NGS-data/Claire_Parage_Terpfactory

→ Check des data

Eurofins a appliqué une "digital normalization" : se renseigner + la dessus.

Installation des outils du pipeline Trinotate / Transdecoder / RSEM avec l’aide de Timothée

Pipeline

[https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification)

[https://trinotate.github.io/](https://trinotate.github.io/)

Running Trinity

[https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)

Evaluation Assemblies

[https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies)

Pour l’analyse des sRNA:

1. Mapping sur la séquence du virus TuMV		OK

2. Unmapped.bam → unmapped.fastq → unmapped.fasta OK au cas où

3. Changé les U en T dans le fichier des MIR → mature.fa, sinon pas de mapping du tout

4. On réduit le set avec** cd-hit-est **

Mise en place de script pour DESeq2 et EdgeR.

Pour le projet de Claire, probablement utiliser RSEM et EBSeq qui permettent de calculer l’expression en prenant en compte tous les isoformes. (nouvel algo avec introduction des pseudo-compte: meilleure estimation des isoformes apparement!)

**Jeudi 6 avril 2016**

Avec l’aide Timothée, mise en place des banques de données et de tous les logiciels pour faire tourner le pipeline de Trinotate

Création blast database avec les assemblages *de novo* de Meja_NT. On peut les query avec "sequenceserver" en interface graphique.

**Mardi 19 avril**

Claire PROJECT:

# Installation de wget sous macOS

brew install wget

# Download des P450 database

wget [http://www.p450.kvl.dk/TFAfiles/AtP450PROT.tfa](http://www.p450.kvl.dk/TFAfiles/AtP450PROT.tfa)

**Message de Stéfanie:**

Pour ce qui est du heatmap,

Pour nicotiana j'avais mappé les smalls contre les séquences de miRNA connus avec 0 mismatch pour distinguer ce type de small (sinon j'avais juste beaucoup trop de données et pas d'info d'annotations)

J'ai calculé pour chaque miRNA son expression en RPM pour chaque échantillon et j'ai calculé des fold change en faisant simplement les rapports de ces expressions entre échantillons.

Les miRNA identifiés ont ensuite étés utilisés pour faire de la prédiction de cibles sur nicotiana à l'aide de psRNAtarget (en ligne). Et j'ai ensuite récupéré pour chaque cible ses données d'expression différentielle (fold change) dans les fichiers de résultats de DESeq issus de l'analyse RNAseq.

Le heatmap consiste ensuite à dessiner côte à côte les fold change des miRs et de leurs cibles correspondantes pour voir s'il y a des liens entre expression des miRs et de leurs cibles.

Je t'ai mis en pièce jointe une des heatmaps que j'ai donné. Il concerne la comparaison des fold change des miRs(m) et des cibles correspondantes(t) entre les différents échantillons.

Je te met aussi le code R utilisé pour générer les heatmaps avec ou sans clustering.

**Mercredi 20 avril**

Réception des données de RNAseq de *Nicotiana benthamiana*.

Installation de la commande parallel sur slurm:

$ yum install parallel -y

Link des fichiers RNASeq Paired-end

[https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbMOCK_1.fastq.gz](https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbMOCK_1.fastq.gz)

[https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbMOCK_2.fastq.gz](https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbMOCK_2.fastq.gz)

[https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV-MUT_1.fastq.gz](https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV-MUT_1.fastq.gz)

[https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV-MUT_2.fastq.gz](https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV-MUT_2.fastq.gz)

[https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV_1.fastq.gz](https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV_1.fastq.gz)

[https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV_2.fastq.gz](https://s3-eu-west-1.amazonaws.com/ascidea-client-ec2/20160405_RNAseq.plant_Donaire/NbTMV_2.fastq.gz)

# Pour download

$ cat raw_reads_links.txt | parallel --gnu "wget {}"

Ou bien

$ wget -i raw_reads_links.txt

Les données sont à mapper sur le même génome de réference que Stéfanie avait utilisé.

A download ici:

[https://solgenomics.net/organism/Nicotiana_benthamiana/genome](https://solgenomics.net/organism/Nicotiana_benthamiana/genome)

wget -r ftp://ftp.solgenomics.net/genomes/Nicotiana_benthamiana/annotation/Niben101/

Sinon il y’a aussi des génomes de référence ici:

[http://sefapps02.qut.edu.au/benWeb/subpages/downloads.php](http://sefapps02.qut.edu.au/benWeb/subpages/downloads.php)

TODO: Check with Stéf & Livia

**Jeudi 21 avril**

Ecriture d’un rapport des analyses réalisées/en cours pour Manfred

Expression différentielle DESeq2 avec les outputs de RSEM sur les données de Claire

Les gènes candidats de Claire sont tous down régulés quand traité avec hormone.

Utilisation des comptages par gènes et non pas pas isoformes pour l’instant.

J’ai l’impression que l’assemblage de novo d’Eurofins est foireux...

TODO: Faire un récap + check fiabilité RSEM

Un des problèmes de l’assemblage de novo de Eurofins, c’est le nombre d'isoform assez impressionant pour certain locus. Va falloir filtrer à mort...

Ou refaire soit-même les assemblages… 

**Vendredi 22 avril**

Test de miRExpress: le soft à l’air plutôt bien pour décrire les données. 

Faire attention au fichier des MIR précurseurs dans **MIRBase_v21**: il y a des codons ambigües pour certaines espèces.

Il faut également convertir les Uracile en Thymine (**U → T**) pour les aligneurs tels que bowtie

<table>
  <tr>
    <td>A</td>
    <td>Adenine</td>
    <td>T</td>
  </tr>
  <tr>
    <td>G</td>
    <td>Guanine</td>
    <td>C</td>
  </tr>
  <tr>
    <td>C</td>
    <td>Cytosine</td>
    <td>G</td>
  </tr>
  <tr>
    <td>T</td>
    <td>Thymine</td>
    <td>A</td>
  </tr>
  <tr>
    <td>Y</td>
    <td>Pyrimidine (C or T)</td>
    <td>R</td>
  </tr>
  <tr>
    <td>R</td>
    <td>Purine (A or G)</td>
    <td>Y</td>
  </tr>
  <tr>
    <td>W</td>
    <td>weak (A or T)</td>
    <td>W</td>
  </tr>
  <tr>
    <td>S</td>
    <td>strong (G or C)</td>
    <td>S</td>
  </tr>
  <tr>
    <td>K</td>
    <td>keto (T or G)</td>
    <td>M</td>
  </tr>
  <tr>
    <td>M</td>
    <td>amino (C or A)</td>
    <td>K</td>
  </tr>
  <tr>
    <td>D</td>
    <td>A, G, T (not C)</td>
    <td>H</td>
  </tr>
  <tr>
    <td>V</td>
    <td>A, C, G (not T)</td>
    <td>B</td>
  </tr>
  <tr>
    <td>H</td>
    <td>A, C, T (not G)</td>
    <td>D</td>
  </tr>
  <tr>
    <td>B</td>
    <td>C, G, T (not A)</td>
    <td>V</td>
  </tr>
  <tr>
    <td>X/N</td>
    <td>any base</td>
    <td>X/N</td>
  </tr>
</table>


Utilisation de miRExpress:

> Raw_data_parse -i TantoT_1.clean.fastq -o TantoT.miRExpress.unique.count.txt

**# Prendre l’output de Raw_data_parse**

> statistics_reads -i TantoT.miRExpress.unique.count.txt -o stats.fasta.test.txt

**Lundi 25 avril**

Pour DrakarM, qui me pose pose problème avec RSEM, voici les résultats du comptage des séquence unique de small RNA:

EDIT: Ne pas utiliser RSEM pour du small !

Les séquences les + abondantes:

2017478	GGGATTGTAGTTCAATTGGTCAGAGCACCGCCC

1925397	GGGATTGTAGTTCAATTGGTCAGAGCACCGCCCT

1544690	GGGATTGTAGTTCAATTGGTCAGAGCACCGCCCC

1101703	AGGGATATAACTCAGCGGTAGAGTGTCACCT

819490	AGTTACTAATTCATGATCTGGC

797755	TGAGGCATCCTAACAGACCGGTAGACTTGAAC

601289	AGGCATCCTAACAGACCGGTAGACTTGAAC

552700	CTAACAGACCGGTAGACTTGAAC

532350	TCCGATGTCGTCCAGCGGTTAGGATATCTGGC

507782	GGGATTGTAGTTCAATTGGTCAGAGCACCGCC

455080	GGGATTGTAGTTCAATTGGTCAGAGCACC

450968	GAGGCATCCTAACAGACCGGTAGACTTGAAC

396933	CAGCTGAGGCATCCTAACAGACCGGTAGACTTGAAC

389004	GGGGATATAGCTCAGTTGGTAGAGCTCCGCT # pas miR

367124	GGTGGCTGTAGTTTAGTGGTGAGAATTCCACGTT # pas miR

348998	CTGAGGCATCCTAACAGACCGGTAGACTTGAAC

291955	GCATCCTAACAGACCGGTAGACTTGAAC

281390	TCCTAACAGACCGGTAGACTTGAAC # pas miR

260117	CATCCTAACAGACCGGTAGACTTGAAC 

258240	GCGTCTGTAGTCCAACGGTTAGGATAATTGCC

236778	ATCCTAACAGACCGGTAGACTTGAAC

Tout ce qui est surligné est du génome chloroplastique… Faudra voir l’importance

Les smalls les plus abondant n’ont pas été trouvé sur mirBASE, du coup j’ai blasté les séquences pour voir. On tombe sur le génome chloroplastique avec 100% d’identité.

![image alt text](image_5.png)

Pour TantoM, pareil… pas mal de chloro, mais également présence de rRNA en bonus…

Les outputs de miRExpress sont plutôt pas mal.

Voir ma présentation au labmeeting de Manfred:

![image alt text](image_6.png)

Mise en place de scripts pour lancer:

Tophat2	OK

Cutadapt	OK

miRExpress	OK		"Attention à bien donner un fichier FASTQ en input !!!"

**For MIR target prediction:**

* **Diana Micro-T**: excellent server for miRNA target prediction in silico, with a friendly user interface. The user can predict targets using miRNAs or putative target gene sequences as input.

* **MirPath**: data mining server for the analysis of pathways controlled by a subset of miRNAs. It is the first server to combine biochemical data with miRNA regulation and the analysis of both factors together constitutes a huge source of information about the implication of miRNAs in biological processes. Very nice interface and easy to interact with.

* **Microinspector**: server for the analysis of potential target sites for miRNAs. The user can input the sequence of the 3'-UTR of the gene of interest, and the server will show the potential miRNA binding sites using the most recent MiRbase version. Several parameters for the target binding can be customized by the user.

* **MaMi**: classical server for the analysis of hybridization energies and structures between miRNAs and target sequences. Chemical parameters can be modified by the user.

* **miRecords**: integrative resource for the target analysis in silico. It combines results from other applications and servers together with bibliographical information about target prediction and validation.

* **mirDIP**: integrates twelve microRNA prediction datasets from six microRNA prediction databases, allowing users to customize their microRNA target searches. Combining microRNA predictions allows users to obtain more robust target predictions, giving you more confidence in your microRNA targets.

* **miRGator**: is a web tool for functional interpretation of miRNAs, integrating functional analysis, expression profile and target prediction to infer physiological roles of miRNAs.

* **Magia**: integrated gene analysis tool for the cross-determinations of miRNA, targets and gene expression. The server allows the analysis of miRNA and genes expression profiles by adopting different statistical measures of profiles relatedness and algorithms for expression profiles combination.

* **Targetscan**: classical software for miRNA target prediction. The server allow to perform searches by miRNA or target gene.

* **miRTar**: MicroRNA Target prediction (miRTar) is a tool that enables biologists easily to identify the biological functionsregulatory relationships between a group of known/putative miRNAs and protein coding genes. It also provides perspective of information on the miRNA targets on alternatively spliced transcripts.

* **miRTrail**: allows you to easily analyse for potential relationships between a set of miRNAs and a set of mRNAs. This enables you to assess possible important implications of the miRNAs on the given disease.

**Mardi 26 avril**

We should discussed about the sRNA distribution and also about the presence of abundant chloroplastic genome in mock samples.

 

**Mercredi 27 avril**

Réunion avec Manfred.

Discussion de l’analyse des smallRNA et de la suite des analyses Brassica napus, ainsi que l’analyse des données de Nicotiana benthamiana.

Download de toute les infos de N.bentha sur le ftp de solgenomic.

$ wget -r ftp://[ftp.solgenomics.net/genomes/Nicotiana_benthamiana/assemblies/](http://ftp.solgenomics.net/genomes/Nicotiana_benthamiana/assemblies/)

**Jeudi 28 avril**

Cutadapt sur les RNAseq Nicotiana benthamiana

Tophat2 mapping en cours 2 mismatches 

miRExpress script au point pour d’autres organismes.

Check de seqbuster/mirAligner

Test de srnaworkbench

![image alt text](image_7.png)

Les paramètres pour psRNAtarget, on a réduit l’expect à 2 pour avoir moins de faux positif

Par contre va falloir trouver une autre méthode pour identifier les target à partir de notre génome à nous.

**Lundi 02 mai**

**Quast**: tool pour évaluer/comparer des assemblages de novo. 

[http://bioinf.spbau.ru/en/quast](http://bioinf.spbau.ru/en/quast)

**Mardi**** 17 mai**

Dear David,

on Monday I start a two weeks vacation; thus, I will be away when you come back. Nevertheless I will read e-mails and also will have time to work in between. 

Since we will have a GAMAVIR network meeting on June2/3 it would be good to advance the analysis of the DE analysis and the correlation with small RNAs as much as possible before the meeting. 

Although we have no replicates yet, we should describe existing data as good as possible. I try below to summarize my current vision of data analysis. However, I do not know what can be done with available tools and time nor whether my vision is the optimal one. Nevertheless, it would be great if you and Khalid could collaborate on this, try a deep analysis of the data, and inform me about your progress by e-mail while I am away. 

1. Genes differentially expressed between mock and infected for both Tanto and Drakkar, and to identify genes commonly regulated. It will be important to use normalized data and to identify a suitable cutoff for the fold change level. Venn diagrams for different cutoff levels would be useful. Moreover, genes with changes in expression (for a specific cutoff level or different cutoff levels) should be annotated and categorised (e.g. COG database, Blast2GO, MapMan). It may also be useful to perform a GO functional enrichment analysis and KEGG pathway analysis. It would be great to get a clue about specific pathways affected by infection. However, I do not know which tools we have in this direction. 

J’ai rajouté les différents informations d’annotation dans le fichier pour Khalid.

J’ai également récup la liste des GoTerm, et converti les ID en texte.

 

2. miRNA levels. It could be useful to make a heatmap (and table with numbers) for the expression of the B napus miRNAs in the 4 different samples. The heat map may allow us to visualize obvious changes upon infection and between the two cultivars for specific miRNAs.

Récup les expressions des différents MIR.

Pour l’instant j’ai identifié les miR de Brassica napus, oleracea, rapa.

Balancer les données dans R pour tout visualiser

3. miRNAs and targets. We spoke about a heatmap correlating the levels of miRNAs with the expression level of known and predicted targets. I have no clear imagination yet how the heatmap could be sorted to visualize important changes and correlations. I assume the data will tell us how to sort them best.  

Les targets ont été identifié avec psRNAtarget en utilisant l’assemblage de novo du génoscope.

Va falloir trouver une nouvelle méthode une fois que les données ADN seront dispo. On utilisera le génome de chaque sample avec les SNPs inférés.

4. the sequenced sRNAs do not belong only to miRNAs. It has been shown that plant produce also other types of sRNAs (e.g. from endogenous IR loci, perhaps annotated in B napus) and that viruses induce the production of virus-induced siRNAs derived from coding genes (vasiRNAs; these are no miRNAs, derived from coding genes, and absent in mock). Thus, if possible, sRNAs should be categorized according to their types and origins. The different sRNA categories should be adressed as under point 2 and 3.

Pour identifier les sRNA qui ne sont pas des MIR

5.  DE analysis starting from degradome data. The degradome data should indicate the sites on cleaved mRNAs corresponding to the presence of corresponding small RNAs in the sRNAseq data. There are different scores (usually 0-4) and I do not exactly know what they mean. For me it is important that the cleaved 5’ end and corresponding sRNA occur with a considerable frequency and that the cleavage site corresponds to the 10-11 position of the sRNA. The number of cleaved targets occurring in the different samples could be shown in venn diagrams (one could show different venn diagrams for different scores).   

The cleaving sRNAs and their targets (with expression level) should be shown in a table/heatmap sorted according to type of sRNA (viral siRNAs, miRNAs, and endogenous siRNAs (i.e. virus-induced siRNAs derived from coding genes, vasiRNAs)) and their frequency/level.

This is a lot of work and we will have to repeat it two times (a. when we get the genome sequence info; b. when we get data of replicates). However it is urgent to develop the analysis pipeline and identify the general trends using this first set of data, so that we can discuss them during the GAMAVIR meeting and proceed efficiently with the new data expected later this year. There will also be soon another ANR review of our project in Paris (I assume in September) and we need data to present (we are already greatly delayed; we should have done already functional analyses of candidates at that time).

Thanks for your help,

Manfred 

Denovo assembly

[http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/)

### Transcriptome assembly

Data used for assembly corresponded to the ~145 million bp of sequence reads generated previously [[3](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR3)], and 293 million bp of new data from 11 Illumina runs covering five tissue-specific libraries. Prior to assembly, the four datasets (thermal-based, head, antennae, and accessory gland) were concatenated, and read abundance was normalized to 50X coverage using the *in silico* normalization tool in Trinity to improve assembly time and minimize memory requirements. Filtering and normalization reduced the dataset to 15 Gb, comprising approximately 32 million normalized read pairs, which were then assembled using default parameters in Trinity (r2014_07-17). Transcript expression levels were estimated with RSEM [[6](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR6)] and open reading frames (ORFs) were predicted using Transdecoder [[7](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR7)]. Hmmer3 was used to identify additional ORFs matching Pfam-A domains. Following transcriptome assembly, reads were filtered, sorted, and prepared for NCBI transcriptome shotgun assembly (TSA) submission as previously described [[8](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR8)].

### Annotation

Functional annotation was performed at the peptide level using a custom pipeline [[8](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR8)] that defines protein products and assigns transcript names. Predicted proteins/peptides were analyzed using InterProScan5, which searched all available databases including Gene Ontology (GO) [[9](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR9)]. BLASTp analysis of the resulting proteins was performed with the UniProt Swiss Prot database (downloaded 11 February 2015). Annie [[10](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR10)], a program that cross-references SwissProt BLAST and InterProScan5 results to extract qualified gene names and products, was used to generate the transcript annotation file. The resulting .gff3 and .tbl files were further annotated with functional descriptors in Transvestigator [[8](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730634/#CR8)].

[http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088462](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088462)

### Sequence Data Analysis and Assembly

The raw data generated by Illumina sequencing were converted from the BCL format to qSeq using Off-line Basecaller, v.1.9.4 (OLB) software. The qSeq files were transformed in FastQ files, which contain sequences that are 72 bp in length, using a custom script. Low-quality sequences were removed; these sequences included reads with ambiguous bases, reads with less than 70 bases, and reads with a Phred quality score Q≤20 using the NGS QC toolkit [[31]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088462#pone.0088462-Patel1). All reads were deposited in the National Center for Biotechnology Information (NCBI) database and can be found under accession number SRA073690.

All datasets were combined, and the sequenced reads were assembled using Trinity ([http://trinityrnaseq.sourceforge.net/](http://trinityrnaseq.sourceforge.net/)), which is a program developed specifically for *de novo*transcriptome assembly from short-read RNA-Seq data that recovers transcript isoforms efficiently and sensitively using the de Bruijn graph algorithm [[32]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088462#pone.0088462-Grabherr1). The optimal assembly results were chosen according to an evaluation of the assembly encompassing the total number of contigs, the distribution of contig lengths, the N50 statistic and the average coverage. The assembled transcripts were based on the main isoform of each transcript, and only contigs with lengths of greater than 300 bp were included in the downstream analysis.

To identify the genotypic contribution to each transcript, reads from each library were mapped against the assembly generated from all libraries using the bowtie aligner [[33]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088462#pone.0088462-Langmead1). The BAM files generated by bowtie were then used to estimate the transcript-level abundance for each library using the RSEM (RNA-Seq by Expectation Maximization) software [[34]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088462#pone.0088462-Li2).

**GO term Quick Tips:**

[http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003343](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003343)

**GO term visualization:**

ftp://[ftp.geneontology.org/pub/go/www/GO.tools_by_type.visualization.shtml](http://ftp.geneontology.org/pub/go/www/GO.tools_by_type.visualization.shtml)

Pour identifier les targets des mir, j’ai utiliser psRNATarget avec les mir de bol, bra et bna.

J’ai pris un expect de 2.0 (au lieu de 3.0) pour réduire le nombre de target. Du coup pour certain MIR il n’y aucun target d’identifier. Peut-être quand même utiliser.

Comparaison des assemblages de novo:

Detonate	[http://deweylab.biostat.wisc.edu/detonate/](http://deweylab.biostat.wisc.edu/detonate/) 

Transrate	[http://hibberdlab.com/transrate/getting_started.html](http://hibberdlab.com/transrate/getting_started.html)

Pipeline pour analyse des assemblages denovo.

Description utilisation de Trinity:

# A basic recommendation is to have ~1G of RAM per ~1M pairs of Illumina reads.

# Simpler transcriptomes (lower eukaryotes) require less memory than more

# complex transcriptomes such as from vertebrates.

# If you are able to run the entire Trinity process on a single high-memory

# multi-core server, indicate the number of butterfly processes to run in

# parallel by the --CPU parameter.

# Our experience is that the entire process can require ~1/2 hour to one hour

# per million pairs of reads in the current implementation.

# Trinity works in a K-mer (K=25) space

# We used the latest version of trinity (v2.2.0)

## 4 stages of Trinity:

# 0. Jellyfish  --> Extracts and counts K-mers (K=25) from reads

# 1. Inchworm   --> Assembles initial contigs by "greedily" extending sequences with most abundant K-mers

# 2. Chrysalis  --> Clusters overlapping Inchworm contigs, builds deBruijn graphs for each cluster, partitions reads between clusters

# 3. Butterfly  --> resolves alternatively spliced and paralogous transcripts independently for each cluster (in parallel)

## Workshop / Tuto

# http://cbsu.tc.cornell.edu/lab/doc/Trinity_workshop_Part1.pdf

# https://github.com/trinityrnaseq/KrumlovTrinityWorkshopJan2016/wiki

# http://dendrome.ucdavis.edu/ftp/Tutorials/PDFs/Trinity_Tutorial_RNA-Seq.pdf

You can parallelized srun job into an sbatch job. 

To see if your srun job are launched, use $ squeue -s

**BLAST INFORMATION:**

Lors d’un blast, on peut ajouter les vrais annotations du* subject *

*Exemple:**  *

gi|1040064781|ref|NR_137327.1| **devient**

Solanum lycopersicum tomato | Solanum lycopersicum 17S ribosomal RNA gene (LOC107882131), ribosomal RNA

Pour ce faire, il faut d’abord download la database **taxdb** avec le script **update_blast.pl**, puis rajouter le path du dossier de la database dans le bash_profile:

export BLASTDB=/le/path/du/dossier/taxdb/:$BLASTDB

Ensuite on peut modifier les outputs du blast, comme sur cette exemple:

blastn -task blastn-short -max_target_seqs 5 -db nr -query NbMOCK_1.clean.fastq.uniq.fasta -remote -out NbMOCK_1.clean.fastq.uniq.fasta.out2 -outfmt "6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle"

Les données de séquençage ADN de* Brassica napus *ont été téléchargé depuis le serveur ftp du BGI ([http://cdts-wh.genomics.cn/](http://cdts-wh.genomics.cn/))

Sauvegardé sur **/mnt/data/bioimage**

### Small RNA processing

Raw sequence reads quality check was performed using FastQC (version 0.11.3). For adaptor trimming cutadapt (version 1.2.1) was used with the following parameters: minimal sequence length after adaptor removal was 16, maximal sequence length was 28 [[65](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR65)]. We used UEA Small RNA Workbench (version 2.5.0) for all miRNA processing tasks [[66](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR66)]. MiRProf calculated the expression of known miRNAs. This program uses Mirbase (version 20). We allowed 2 mismatches during the analysis. We grouped together mismatches, variants and miRNAs from different organisms and worked only with families. MirCat identified candidant microRNAs. We used the default plant specific parameters. PAREsnip was used to process degradome sequences for finding microRNA targets [[40](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR40)] (Additional file [11](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#MOESM11)). During this analysis, we discarded results from category 4, the maximal mismatch number was 4, mismatches at position 10–11 were not allowed and more than two mismatches at adjacent positions also not allowed. Every other parameter was the default. We have created a non-redundant scaffold set from the available *N. benthamiana* genome sequencing projects using Blast. We used this set as genome in all programs where genome sequence was a mandatory requirement. Results from these programs were analysed further with custom made Perl and Python scripts. For statistical analysis and visualization R was used [[67](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR67)].

### Transcript processing and annotation

To generate the corresponding mRNA transcriptome we have sequenced Illumina paired end TruSeq libraries from the same *N. benthamiana* leaf, stem, germ, root, flower and seed samples as at the small RNAs. We pooled together 315 million 2x100 bp reads and used Trinity *de novo* assembly program (r2013-02-16) with default parameters [[68](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR68)].

To find biological function for de-novo transcripts, Trinotate was used (version: 20130225) according to the manual [[68](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR68)]. All related data were pushed into an SQLight database for further processing.

To perform the Gene Ontology (GO) analysis of miRNA targets the agriGO web-based tool and database version 1.2 was used [[69](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2209-6#CR69)]. The analysis setting was the following: Singular Enrichment Analysis (SEA) with the *Arabidopsis thaliana* TAIR9 reference gene model. The graphical results showing the molecular functions.

RNAseq analysis benthamiana replicate 1 (1 zone)

**Common up-regulated (>2FC) gene in infected vs mock:**

**Common down-regulated (<-2FC) gene in infected vs mock:**

[https://www.araport.org/downloads/Araport11_Release_201606/annotation](https://www.araport.org/downloads/Araport11_Release_201606/annotation)

Araport va continuer le support sur Arabidopsis et remplacer TAIR10 du coup.

# Graphical soft showing the sRNA repartition on virus after mapping: **LayerCake**

[http://graphics.cs.wisc.edu/WP/layercake/how-to-use-layercake/](http://graphics.cs.wisc.edu/WP/layercake/how-to-use-layercake/)

Dear all,

We came a long way in collecting data and data analysis. However, the accumulating data and associated questions and decisions to take are also becoming overwhelming and this puts us into a situation in which we need more organization and transparency, and a platform for discussion. 

Thus, I discussed with Stéfanie to have a seminar-style meeting with presentations about the "TMV-benthamiana" and "napus-TuMV (GAMAVIR)" projects, in which the status quo of the analyses and the further steps should be presented and discussed. Presentations may show preliminary tables for discussion but should also start to focus on the presentation of final (and preliminary) results in the form of fully annotated summarizing slides (clear to us but also to other audiences). The goal of the meeting shall be to 

(1) determine what has been achieved - distinguishing between final data and preliminary data -,  

(2) present/set the goals for further bioinformatic and wetlab analysis, reports and publication, and 

(3) discuss the formats of final data output and visualization.

It is important to decide which bioinformatic data are final, which ones to keep for further analysis, and which ones to put aside. 

Stéfanie and David should present their workflow and data outputs and Khalid and Nicolas should present what they did, do and will do with these data outputs. With the goals in mind we should discuss together the next steps and who will do what, and thereby optimize the teamwork. This seminar should become the first of a series of seminars to discuss the project progress in regular terms. This interaction format should be of mutual benefit and enhance our collaboration. I hope that most of the bioinformatics analysis can be completed by the end of the year so that we can then focus on experiments and publication. I also need to present GAMAVIR results during Plant-KBBE auditing (likely) in December, which sets an important deadline. 

I assume that our (first) meeting in this format will take a long time (4-6 hours?). With Stéfanie I propose the 12th of October and I already reserved the seminar room 243 for 11:00-18:00 for us. If this date is not possible, please let me know. 

All the best,

Manfred

Cutadapt option for sRNA cleaned →

cutadapt_options="cutadapt -a TGGAATTCTCGGG --trim-n --minimum-length 18 --maximum-length 30 -q 30 --discard-untrimmed -o $cleanreads_dir/${read_id}.clean.fastq.gz $rawreads_dir/${read_id}.fastq.gz"

* Re-clean les librairies en enlevant les 18-30 !

