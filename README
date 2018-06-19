# Installation and usage documentation of FUNN-MG pipeline version 1.12
# Last edited: 26/01/2017
# Author: Leandro Corrêa

CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Software requirements
    * Platform
    * Third-party software requirements
 * Installation
 * Configuration
 * Input files
 * How to run
 * Output files
 * For docker users
    * Execution in Linux Environment
    * Execution in OS X Environment
 * Troubleshooting
 * FAQ
 * Maintainers

INTRODUCTION
------------

FUNN-MG is a tool for functional and visual analysis of bi-partite networks (genes and metabolic pathways) calculated from metagenomics data. The pipeline is written in R version 3.3.1.

SOFTWARE REQUIREMENTS
----------------------

----------------------
PLATFORM
----------------------

    * Linux, (32 bit or 64 bit) or MAC OS X Mavericks (64 bit recommended)

    * R version >= 3.2.1.

    * Java JRE/JDK version 6 or higher (Oracle or OpenJDK) for RedeR package.

    * libcurll. For linux users: sudo apt-get install libcurl4-gnutls-dev.

    * libxml. For linux users: sudo apt-get install libxml2-dev.
	
    * libssl-dev. For linux users: sudo apt-get install libssl-dev.

    * libsasl2-dev. For linux users: sudo apt-get install  libsasl2-dev.

    * pymongo. For linux users: sudo python -m pip install pymongo

    * pandas. For linux users: sudo python -m pip install pandas

----------------------
THIRD-PARTY SOFTWARE REQUIREMENTS
----------------------

All R packages dependencies can be found in folder packages, but they are automatically installed with the FUNN-MG.

    * gtable version >= 0.2.0

    * qvalue version >= 2.4.2

    * RCurl version >= 1.95-4.8

    * RedeR version >= 1.20.0

    * reshape2 version >= 1.4.1

    * scales version >= 0.4.0

    * XML version >= 3.98-1.4


Note: Internet connection are required to install packages and run queries on KEGG database.


INSTALLATION
------------

   1.) Install Java, libcurl and libxml dependences as defined under the section SOFTWARE REQUIREMENTS -> PLATFORM

   2.) Download the FUNN-MG pipeline

       2.2) Uncompress and untar the archive file:

            tar -vzxf funn-mg-v.*.tar.gz

   Warning: The first execution of the FUNN-MG may take a few minutes because of the installation of R packages dependencies.


CONFIGURATION
------------

INPUT FILES
------------

FUNN-MG pipeline takes as input a list of ko groups in csv format, containing the identifiers of the sequences and homology scores.

    +------------------------------------------------------------------------+-----------------------------------------------------------------------+
    |                                                            Input file model - FUNN-MG                                                          |
    +------------------------------------------------------+-----------------+-----------------------------------------------------------------------+
    |        ID        |      KO        |        GENE      |       SCORE     |                               SEQUENCE                                |
    +------------------+----------------+------------------+-----------------+-----------------------------------------------------------------------+
    |   id00341        |    K06168      |    bth:BT_3195   |      199.519    |  dbj|AP014924.1| Limnochorda pilosa DNA, complete genome, strain:...  |
    +------------------+----------------+------------------+-----------------+-----------------------------------------------------------------------+
    |   id00441        |    K00266      |    bth:BT_4310   |      1336.24    |  dbj|AP012211.1| Eggerthella sp. YY7918 DNA, complete genome          |
    +------------------+----------------+------------------+-----------------+-----------------------------------------------------------------------+
    |   id00634        |    K02986      |    bth:BT_2702   |      109.768    |  gb|CP008876.1| Terribacillus aidingensis strain MP602, complete ...  |
    +------------------+----------------+------------------+-----------------+-----------------------------------------------------------------------+
    |   id00635        |    K02986      |    bth:BT_2702   |      72.4034    |  ref|XM_008022009.1| Setosphaeria turcica Et28A hypothetical prot...  |
    +------------------+----------------+------------------+-----------------+-----------------------------------------------------------------------+
    |   id00903        |    K02945      |    bth:BT_4345   |      1105.89    |  gb|CP003470.1| Rhodanobacter denitrificans strain 2APBS1, comple...  |
    +------------------+----------------+------------------+-----------------+-----------------------------------------------------------------------+


Note: SCORE, ID, GENE and SEQUENCE are optional parameters.



HOW TO RUN
------------

----------------------
COMMAND LINE ARGUMENTS
----------------------

   1.) In the terminal, go to the src folder:

       cd ~/funn-mg-v*/src

   2.) Run example: Rscript funn-mg.R -i <[INPUT-FILE-PATH]> -o <[OUTPUT-DIRECTORY-PATH]> -p <[NUMBER]> -t <[CHARACTER]> -f <NOT VALUE> -g <NOT VALUE> -d <NOT VALUE>

        -i
            (required) Input file path (relative or absolute).

        -o
            (required) Defines an existing absolute path where an output subdirectory is written to.

        -p
            (optional) Set the cutting range in the confidence test of statistical coverage [0..1]. Defaut = 0.05.

        -t
            (optional) View KEGG pathways by hierarchy: [c] Classe; [s] Subclass; [p] ID pathwway. Default = [p].

        -f
            (optional) Fast execution. It only works if the pipeline runs at least once.

        -g
            (optional) Consider global and overview maps in the analyse.

        -d            
	    (optional) Disable network view.


Note: For more information run Rscript funn-mg.R --help


OUTPUT FILES
------------

1.) FUNN-MG displays a interactive network of elements identified and selected in the sample. That view is part of RedeR package and allows among other functions: zoom, pan, query (ctrl + F), addition and deletion of network elements, grouping...


2.) FUNN-MG also includes a .csv file with a summary of the analysis along with topological network metrics identified for each pathway annotated in the sample:

   * pathways: id that identify the metabolic pathyway in the KEGG database.
     --------

   * class: class of metabolic pathway by the KEGG hierarchy.
     ------

   * subclass: subclass of metabolic pathway by the KEGG hierarchy.
     --------

   * Function: biochemical functions corresponding to each pathway.
     ---------

   * genes_relation_sample: number of genes found in the sample involved in the pathway.
     ---------------------

   * genes_relation_noted: number of genes noted in the KEGG database involved in the pathway.
     --------------------

   * coverage: hypergeometric test applied to the pathways coverage analysis. While lower results, greater the chances the pathways have not been identified at random.
     --------

   * p.value: result of Fisher exact test for confidence index of coverage analysis.
     -------

   * q.value: result of FDR test for p.value index.
     -------

   * betweenness_centrality: It is an indicator of a pathway's centrality in a network.
     ----------------------

   * neighborhood_connectivity: It is defined as the average connectivity of all genes neighbors of the pathway.
     -------------------------

   * funn_pathways: metabolic pathways selected by the FUNN-MG.
     -------------


3.) If there is an information about the species related to each KO group in the input table, FUNN-MG also includes a .txt file containg the species related to each path mapped in the network, in descending order according to network centrality scores.



FOR DOCKER USERS
----------------

-----------------------------
EXECUTION IN LINUX ENVIROMENT
-----------------------------

1.) Install docker: https://docs.docker.com/engine/installation/linux/


2.) Push hscleandro/funn-mg docker image

        ~> sudo docker push hscleandro/funn-mg


3.) Access the docker image (hscleandro/funn-mg). Example:

        ~> xauth list
           itvdsws058/unix:1  MIT-MAGIC-COOKIE-1  9bb1b83b379efa91d5c3fbc412a9c9ef
           itvdsws058/unix:0  MIT-MAGIC-COOKIE-1  17b88f1c67afb0c0b7a6b58316a90d5c
           ITVDS-WS058/unix:0  MIT-MAGIC-COOKIE-1  58c68d6c6af0b6354d66cd9581c5d9eb

    * Copy the last row of output.

        ~> sudo docker run -t -i --net=host -e DISPLAY -v /tmp/.X11-unix/ hscleandro/funn-mg /bin/bash

    * Inside docker container add the last row copied previously to the xauth command. Example:

        ~> xauth add ITVDS-WS058/unix:0  MIT-MAGIC-COOKIE-1  58c68d6c6af0b6354d66cd9581c5d9eb

4.) Go to the /home/funn-mg*/src diretorie and run.

   * Note: It is still important to install java > 6 (outside the docker) to view the network.


-----------------------------
EXECUTION IN OSX ENVIRONMENT
-----------------------------

1.) Install docker via Homebrew : http://brew.sh/
	 
	 ~> brew install docker
	    brew install docker-machine
	    brew cask install xquartz 


2.) Open XQuartz in terminal mode (it can be on iTerm2) and set the display to your IP (ifconfig on OS X)

	 ~> open -a XQuartz
	    export DISPLAY="172.2.0.74:0.0"
	    ip=$(ifconfig en1 | grep inet | awk '$1=="inet" {print $2}') 
	    xhost + $ip

3.) Run docker
	 ~> docker-machine start default
	    eval $(docker-machine env default)
 
4.) Access the docker image (hscleandro/funn-mg). Example:

	    xauth list

    * Copy the last output.

        ~> docker run -t -i --net=host -e DISPLAY=$ip:0 -v /tmp/.X11-unix:/tmp/.X11-unix hscleandro/funn-mg /bin/bash

    * Inside docker container add the last row copied previously to the xauth command. Example:

        ~> xauth add 172.2.0.74:0  MIT-MAGIC-COOKIE-1  94040ffedd267e5587774173c7673a28

5.) Go to the /home/funn-mg*/src diretorie and run.

    * to run funn analysis and plot the network. Example:

        ~> Rscript funn-mg.R -i ../input/ITV/Lagoas/TRES_IRMAS_BLAST_KAAS.csv -o ../output/ITV/ -t c

    * to plot only the network (it had already did the funn analysis). Example:

        ~> Rscript funn-mg.R -i ../input/ITV/Lagoas/TRES_IRMAS_BLAST_KAAS.csv -o ../output/ITV/ -t c -f

   * Note: It is still important to install java > 6 (outside the docker) to view the network.


Troubleshootinginer
------------

...

FAQ
------------

...

Maintainers
------------

Current maintainers:

* Leandro Corrêa (ITV-DS/UFPA) - hscleandro@gmail.com
* Ronnie Alves   (ITV-DS/UFPA) - alvesrco@gmail.com
nstallation and usage documentation of FUNN-MG pipeline version 1.12
# Last edited: 1/02/2017
# Author: Leandro Corrêa, Ronnie Alves

