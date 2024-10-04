# Zukunfsorientiertes Workflow-Design für reproduzierbare Ergebnisse von ATACseq und RNAseq Analysen

Mithilfe der in diesem Repository enthaltenen Pipelines können die größten Teile der in der Publikation [Variation in Pleiotropic Hub Gene Expression Is Associated with Interspecific Differences in Head Shape and Eye Size in Drosophila](https://doi.org/10.1093/molbev/msaa335) [1] beschriebenen Analysen durchgeführt werden.
Dieses Projekt soll zeigen mit welchen Mitteln ein Analyse-Workflow einer wissenschaftlichen Arbeit zukunftsorientiert gestaltet werden kann, um reproduzierbare Ergebnisse generieren zu können. Weitere Informationen gibt es dazu in der gleichnamigen Arbeit.

## Wie benutze ich die Pipelines?

### Systemanforderungen
Für diese Pipelines wird folgende Hardware als Mindestanforderung empfohlen:

- 8-Core 3 GHz CPU mit Multithreading (besser: 12 Core 4GHz)
- 32 GB RAM
- \>500 GB freier Speicherplatz

### Abhängigkeiten

- [git](https://git-scm.com) [2]
- [Conda](https://docs.anaconda.com/miniconda/) [3]
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) [4]
- [Docker](https://docs.docker.com/engine/install/) [5]

### 1. Vorbereitung

- Installiere alle genannten Abhängigkeiten
- Klone dieses Repository
- Öffne das Repository im Terminal
- Bevor Piplines ausgeführt werden können, müssen die zur Verfügung stehenden Ressourcen des Systems in der Konfigurationsdatei `nextflow.config`angepasst werden. Standardmäßig sind die oben genannten Systemanforderungen eingetragen.

### 2. Download der Daten

- Um die benötigten Daten herunterzuladen, steht eine eigene Pipeline zur Verfügung. Diese kann mittels `nextflow run DownloadData.nf` aufgerufen werden. Die Daten werden in einer vordefinierten Datenstruktur für die RNAseq- und ATACseq-Pipelines gespeichert. Weitere Informationen zu Parametern finden sich hier.

### 3. RNAseq Pipeline

Die RNAseq Pipeline kann vollständig mit dem Befehl `nextflow run RNAseq.nf --align --DEG` ausgeführt werden. Da das Alignment - abhängig der Systemvorraussetzungen - mehrere Stunden in Anspruch nehmen kann, ist es möglich die Differenzielle Genexpressionsanalyse (`--DEG`) unabhängig vom Alignment auszuführen. Vorraussetzung ist, dass das Alignment (`--align`) zuvor einmal ausgeführt wurde. Weitere Informationen zu Parametern finden sich hier.

### 4. ATACseq Pipeline

#### 4.1 Pnr-Motif

Um mithilde der ATACseq Pipeline de novo ein Pnr-Motif zu generieren, wird der Befehl `nextflow run ATACseq.nf --pnrmotif` aufgerufen. Hierfür werden die im Dateipfad `/data/fasta/ChIP-chip` liegenden BED Daten herangezogen und die Ergebnisse in `/data/results/ATACseq/de-novo_Pnr-motif` gespeichert. Aus diesen Resultaten muss das beste Pnr-Motife manuell herausgesucht und im Ordner `/data/pnr-motif` abgelegt werden. Alle weiteren Pnr-Motife müssen in `/data/pnr-motif/optional` abgelegt werden. Weitere Informationen zu Parametern finden sich hier.

#### 4.2 ATACseq Workflow

Ist erfolgreich ein Pnr-Motif gefunden und abgelegt, so kann die restliche Pipeline ausgeführt werden: `nextflow run ATACseq.nf --align --peakcalling --motif`. Wie in der RNAseq Pipeline, kann auch hier die restliche Pipeline getrennt vom Alignment ausgeführt werden, wenn zuvor das Alignment (`--align`) ausgeführt wurde. (Hinweis: `--peakcalling` kann auch ohne `--motif` ausgeführt werden, allerding muss `--motif` zusammen mit `--peakcalling` ausgeführt werden.) Weitere Informationen zu Parametern finden sich hier.

## Methoden

### DownloadData.nf

```shell
nextflow run DownloadData.nf <params>

Params:
--download_file "/path/to/file.csv" (default: "$baseDir/data/fasta/SRA_Download.csv")
                Absoluter Dateipfad der *.csv Datei, die den zu herunterzuladenden SRA Run, den Ausgabenamen und relativen Ausgabepfad zu --data_path im Dateisystem enthält.
--data_path     "/path/to/folder" (default: "$baseDir/data")
                Absoluter Dateipfad zum übergeordneten Ordner, der alle Daten enthalten wird
--rna_ref       "http://link.to/rna_reference.gz" (default: https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-transcript-r6.54.fasta.gz)
                Referenz Transkriptom für die RNAseq Analyse
--atac_ref      "http://link.to/atac_reference.gz" (default: https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz)
                Referenz Genom für die ATACseq Analyse im FASTA Format
--atac_gtf      "http://link.to/atac_gtf_reference.gz" (default: https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz)
                Referenz Genom für die ATACseq Analyse im GTF Format
--chip          ["link.to/bigwig/chip1","link.to/bigwig/chip1"] (default: ["http://furlonglab.embl.de/labData/publications/2012/Junion-et-al-2012_Cell/signal_files/Pnr_4-6h_allIPvsallM_log2ratio_5probes-smoothed_Junion2012.bigwig", "http://furlonglab.embl.de/labData/publications/2012/Junion-et-al-2012_Cell/signal_files/Pnr_6-8h_allIPvsallM_log2ratio_5probes-smoothed_Junion2012.bigwig"])
                Link Liste zu ChIP-Chip Daten
```

### RNAseq.nf

```shell
nextflow run RNAseq.nf <params>

Params:
--align     Wenn gesetzt, dann wird Alignment der Reads durchgeführt
--DEG       Wenn gesetzt, dann wird Differential gene expression durchgeführt

--genome    "/path/to/*.fasta" (default: "$baseDir/data/references/RNAseq/*.fasta")
            Absoluter Dateipfad zu Referenz-Transkriptom
--reads     "/path/to/*.fasta" (default: "$baseDir/data/fasta/RNAseq/**/*.fasta" )
            Absoluter Dateipfad zu Reads
--outdir    "/path/to/results" (default: "$baseDir/data/results/RNAseq")
            Absoluter Dateipfad in dem Resultate gespeichert werden sollen
```

### ATACseq.nf

```shell
nextflow run RNAseq.nf <params>

Params:
--align         Wenn gesetzt, dann wird Alignment der Reads durchgeführt
--peakcalling   Wenn gesetzt, dann wird Peak Calling durchgeführt
--pnrmotif      Wenn gesetzt, dann wird de novo Pnr-Motif generiert
--motif         Wenn gesetzt, dann wird Pnr-Motif in Peaks gesucht

--genome        "/path/to/*.fasta" (default: "$baseDir/data/references/ATACseq/*.fasta")
                Absoluter Dateipfad zu Referenz-Genom im FASTA Format
--gtf           "/path/to/*.gtf" (default: "$baseDir/data/references/ATACseq/*.gtf")
                Absoluter Dateipfad zu Referenz-Genom im GTF Format
--reads         "/path/to/*.fasta" (default: "$baseDir/data/fasta/ATACseq/**/*.fasta" )
                Absoluter Dateipfad zu Reads
--chip          "/path/to/*.bed" (default: "$baseDir/data/fasta/ChIP-chip/*.bed")
                Absoluter Dateipfad zu ChIP-chip Daten
--pnr_motif     "/path/to/motif/*.motif" (default: "$baseDir/data/pnr-motif/*.motif")
--pnr_motif_opt "/path/to/opt/motifs/*.motif" (default: "$baseDir/data/pnr-motif/optional/*.motif")
--outdir        "/path/to/results" (default: "$baseDir/data/results/ATACseq")
                Absoluter Dateipfad in dem Resultate gespeichert werden sollen
```

## Quellenverzeichnis

> [1] Buchberger, et al.**Variation in pleiotropic hub gene expression is associated with interspeciﬁc diﬀerences in head shape and eye si- ze in drosophila.** 
 Molecular Biology and Evolution, 38(5) : 1924–1942, January 2021. doi:[10.1093/molbev/msaa335](https://doi.org/10.1093/molbev/msaa335).

> [2] Scott Chacon und Ben Straub. **Pro Git.** Apress, 2014. doi:[10.1007/978-1-4842-0076-6](http://dx.doi.org/10.1007/978-1-4842-0076-6)

> [3] Anaconda Inc. **Miniconda - Anaconda documentation**. url:https://docs.anaconda.com/miniconda/

> [4] P. Di Tommaso, et al. **Nextﬂow enables reproducible computational workﬂows.** Nature Biotechnology, 35(4) :316–319, April 2017. doi:[10.1038/nbt.3820](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3820.html)

> [5] Docker Inc. **Install Docker Engine.** url:https://docs.docker.com/engine/install/ [zugegriﬀen am 03-10-2024]]