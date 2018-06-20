# FUNN-MG pipeline version 1.13

FUNN-MG is a tool for functional and visual analysis of bi-partite networks (genes and metabolic pathways) calculated from metagenomics data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

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

3.) if you downloaded, uncompress and untar the archive file:
```
tar -vzxf funn-mg-v.*.tar.gz
```
4.) Go to "view" and "funn-mg-v.13" folders and look for the "setup.conf" files (they are two files one in each folder). 

4.1) Enter in the in each setup.conf file and address the "host" and "port" of your mongodb. 

4.1.1) In the field "output" choose a folder that will receive the temporary files of execution that **should not be erased**.

5.) Pending R packages will be installed on first run.

*Warning*: The first execution may be more time-consuming.

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing


## Versioning
FUNN-MG version 1.13
Last edited: 20/06/2018

## Authors

* **Leandro Corrêa** - *Developer, , researchdesigner* - [web page](https://hscleandro.wixsite.com/professional)

* **Ronnie Alves** - *supervisor, research* - [web page](https://sites.google.com/site/alvesrco/)

## License

This project is licensed under the GNU general public licence - see the [LICENSE](funn-mg-v.13/LICENSE.md) file for details

## Acknowledgments

* Instituto Tecnológico Vale Belém/PA [ITV/DS](http://www.vale.com/brasil/PT/initiatives/innovation/itv/Paginas/default.aspx)
* The visualization was inspired by the work of [Dean Atali](https://deanattali.com/).

