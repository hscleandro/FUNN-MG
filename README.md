# FUNN-MG pipeline version 1.13

FUNN-MG is a tool for functional and visual analysis of bi-partite networks (genes and metabolic pathways) calculated from metagenomics data. This release contains the integration of results from [KAAS](https://www.genome.jp/kaas-bin/kaas_main?prog=GHOSTX&way=s), [Kaiju](http://kaiju.binf.ku.dk/), and data derived from the proteomic analysis.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

* Linux, (32 bit or 64 bit) or MAC OS X Mavericks (64 bit recommended)

* R version >= 3.2.1.

* Java JRE/JDK version 6 or higher (Oracle or OpenJDK).

* libcurll. For linux users: sudo apt-get install libcurl4-gnutls-dev.

* libxml. For linux users: sudo apt-get install libxml2-dev.
	
* libssl-dev. For linux users: sudo apt-get install libssl-dev.

* libsasl2-dev. For linux users: sudo apt-get install  libsasl2-dev.

* pymongo. For linux users: sudo python -m pip install pymongo

* pandas. For linux users: sudo python -m pip install pandas

* mongodb >= 3.2. [Guide for Linux users](https://www.digitalocean.com/community/tutorials/how-to-install-mongodb-on-ubuntu-16-04)

* Internet connection.

### Installing

1.) Install the dependences as defined under the section above

2.) Download or clone the FUNN-MG pipeline

3.) If you downloaded, uncompress and untar the archive file:
```
tar -vzxf funn-mg-v.*.tar.gz
```
4.) In "view" and "funn-mg-v.13" folders look for the "setup.conf" files (they are two files one in each folder). 

4.1) In setup.conf files address the "host" and "port" of your mongodb. 

4.1.1) In the field "output" choose a folder that will receive the temporary files of execution. **This folder should not be erased**.

5.) FUNN-MG will install the pending R packages on the first run.

*Warning*: The first execution may be more time-consuming.

### Input files

Funn-mg supports three types of files from different analyzes: KAAS output, Kaiju output, proteomic analysis (amino acids fasta sequence). The ids of each analysis must be the same for each sample.

*Warning*: Only KAAS output data is required for tool execution, the other analyzes are optional.

#### **KAAS output, example:**

```
contig00001_17393_24748_+
contig00001_24760_26394_+	K06160
contig00001_26768_28021_+	K09458
contig00001_28077_28589_+
contig00001_28759_29775_+	K00648
contig00011_4826_5912_-
```
#### **Kaiju output, example:**

```
C	contig00001_17393_24748_+	34078	Bacteria; Cyanobacteria; NA; Nostocales; Scytonemataceae; Scytonema; Scytonema hofmannii; 
C	contig00001_24760_26394_+	1163	Bacteria; Cyanobacteria; NA; Nostocales; Nostocaceae; Anabaena; NA; 
C	contig00001_26768_28021_+	1245922	Bacteria; Cyanobacteria; NA; Nostocales; Scytonemataceae; Scytonema; Scytonema millei; 
C	contig00001_28077_28589_+	28072	Bacteria; Cyanobacteria; NA; Nostocales; Nostocaceae; Nostoc; Nostoc sp. PCC 7524; 
C	contig00001_28759_29775_+	142864	Bacteria; Cyanobacteria; NA; Nostocales; Nostocaceae; Cylindrospermum; Cylindrospermum stagnale; 
U	contig00011_4826_5912_-	0
```
*Warning*: The results of the kaiju tool should contain all taxonomic classifications. For more datails access the [kaiju project](https://github.com/bioinformatics-centre/kaiju/blob/master/README.md).

#### **Hit protein fasta (.faa), example:**

```
>contig00001_17393_24748_+
AQLDLLRHQLSPQEVAGRTRAFIIGGENLVAQTIDFWQEFAMQYEGSIAKTSTQISHKSS
>contig00001_24760_26394_+
STTQHELWLRANQGDAKRLKRPNSVMQIMVYLYSA
>contig00001_26768_28021_+
KNWQLKRVVVTGMGAITPLGNTVTEYWQGLLQGRSGIHPIT
>contig00001_28077_28589_+
KKHVTQNITDPFIGNFNNKIQYFEGVLPEVLSFRIASYQICKNWLKAREGSAFSDEDSHQYNRIVMI
>contig00001_28759_29775_+
LKEIIKLTEEIKTAIQCYHLNKLEIYEKVRAIVVDKLEIEPERVTPTANFSKDLGADSLDTVELVMALEEAFDIEISEQVAKTLLTVQQAIDYISQKVKFAV
>contig00011_4826_5912_-
ATQAAQRAIAMAGLAPKEIDLIILATSTPDDLFGNA
```

## Running the tests

*Warning*:  Before you run the tool, make sure if mongodb is active.

### **Sequence analysis**

1.) In the terminal, go to the src folder:
```
cd ~FUNN-MG/funn-mg-v*/src
```
2.) Download the test file in this [link](https://www.dropbox.com/sh/u31hthujjrnb7m8/AADlfwjfR7dEmacqQOX2qGM0a?dl=0).


3.) Running this command:
```
Rscript funn-mg.R -s EMMA -p LAKES -f '/emma_lakes.kaas' -t '/emma_lakes_kaiju.out' --prot '/emma_lakes_proteomics.fasta'
```

### **View the results**

1.) In the terminal, go to the FUNN-MG folder:

```
cd ~FUNN-MG
```

2.) Running this command:

```
R -e "shiny::runApp('view')"
```

3.) Using your browser for access the address given at the end of the run.


### **Plot network**

1.) In the terminal, go to the src folder:

```
cd ~FUNN-MG/funn-mg-v*/src
```

2.) Running this command:

```
 Rscript funn-mg.R -s EMMA -p LAKES --display --prot
``` 

3.) Running this command:

```
Rscript funn-mg.R -s EMMA -p LAKES --display --type class
``` 

4.) Running this command:

```
Rscript funn-mg.R -s EMMA -p LAKES --display --type subclass
``` 

### Help

For information about the FUNN-MG options run: 
```
Rscript funn-mg.R --help
```

### Video test

You can check this example running in this [video](https://youtu.be/m6aMrSSCvO0).


## Built With

* [R](https://www.r-project.org/) - Main computational language used
* [shiny](https://maven.apache.org/) - The web framework used Dependency Management
* [python](https://www.python.org/) - Support computational language used
* [mongodb](https://www.mongodb.com/) - Database for storage the data processed
* [KEGG](https://www.genome.jp/kegg/pathway.html) - Data base for relationship identification among the omics components


## Contributing


## Versioning

FUNN-MG version 1.13. Last edited: 20/06/2018

## Authors

* **Leandro Corrêa** - *developer, researcher, software architect* - [personal web page](https://hscleandro.wixsite.com/professional)

* **Ronnie Alves** - *supervisor, researcher* - [personal web page](https://sites.google.com/site/alvesrco/)

## License

This project is licensed under the GNU general public licence - see the [LICENSE](https://github.com/hscleandro/FUNN-MG/blob/master/funn-mg-v.13/LICENSE) file for details

## Acknowledgments

* Instituto Tecnológico Vale Belém/PA [ITV/DS](http://www.vale.com/brasil/PT/initiatives/innovation/itv/Paginas/default.aspx)
* The visualization was inspired by the work of [Dean Atali](https://deanattali.com/).

