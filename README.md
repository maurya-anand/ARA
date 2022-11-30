## AMA : An automatic pipeline for exploration of SRA datasets with sequences as a query


### Requirements

-   **Docker**

    -   Please checkout the [Docker installation](https://docs.docker.com/get-docker/) guide.

        *or*

-   **Mamba package manager**

    -   Please checkout the [mamba or micromamba](https://mamba.readthedocs.io/en/latest/installation.html) official installation guide.

    -   We prefer `mamba` over [`conda`](https://docs.conda.io/en/latest/) since it is faster and uses `libsolv` to effectively resolve the dependencies.

    -   `conda` can still be used to install the pipeline using the same commands as described in the installation section.

        > Note: **It is important to include the 'bioconda' channel in addition to the other channels as indicated in the [official manual](https://bioconda.github.io/#usage "Bioconda - Usage")**. Use the following commands in the given order to configure the channels (one-time setup).
        >
        > ``` bash
        > conda config --add channels defaults
        > conda config --add channels bioconda
        > conda config --add channels conda-forge
        > conda config --set channel_priority strict
        > ```

------------------------------------------------------------------------

### Installation

The user can install the pipeline by using either Docker or Mamba using the steps mentioned below.

Firstly, download and extract the zip file into the desired location before starting the setup.

Before starting any analysis with the pipeline, please make sure that the system has enough disk space available for the data you wish to retrieve and process from the SRA repository.

-   **Using Docker**

    ``` bash
    cd AMA/
    docker build -t ama_img .
    ```

*or*

-   **Using Mamba**

    ``` bash
    cd AMA/
    mamba env create --file requirements.yaml
    mamba activate ama_env
    perl setup.pl
    ```

    > *Note: After installation, the virtual environment consumes approximately 1.5 GB of disk space. The installation was tested on "Ubuntu 20.04.4 LTS", "Ubuntu 22.04.1 LTS" and "Fedora 37" using the procedure mentioned above.*

Please be patient because downloading and configuring the tools/modules may take several minutes. The warning messages that appear during the installation of certain Perl modules can be ignored by users.

Refer the [Troubleshooting](#troubleshooting) section in case of any installation related issues.

------------------------------------------------------------------------

### Example usage

-   **Docker**

    `docker run -it ama_img /home/AMA-main/ama.pl --input /home/AMA-main/example/SraRunInfo.csv --sequences /home/AMA-main/example/Arabidopsis_thaliana.TAIR10.ncrna.fa`

-   **Mamba environment**

    `perl ama.pl --input example/SraRunInfo.csv --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa`

To get full usage info: `perl ama.pl --help`

> *Note*: The user can delete the contents of `results/` directory after testing the tool using the example mentioned above.

### Configuration file

The configuration file `conf.txt` is automatically generated during the installation by setup script. It contains certain default parameters as well as the location to the executable binaries of the [tools](#list-of-perl-modules-and-tools-incorporated-in-the-pipeline) incorporated in the pipeline.

The user can modify the default parameters in `conf.txt` and pass it to the pipeline as an input. For example, the user can choose between *blastn* or *bowtie2* by changing the 'execute flag' to either 0 or 1 in the configuration file while leaving the rest of the parameters to default values. By default, both the tools are enabled *ie*. `execute = 1`.

By default, the pipeline uses a pre-built Kraken2 viral genomic database ([release: 9/8/2022](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220908.tar.gz)) from <https://benlangmead.github.io/aws-indexes/k2>. Users can provide their own database by changing the `kraken2_db_path` parameter in the `conf.txt` file.

> *Note:* If the user wishes to use a different installation than Bioconda, the user can manually install the required [tools](#list-of-perl-modules-and-tools-incorporated-in-the-pipeline) and specify the absolute path of the executable binaries in the configuration.

------------------------------------------------------------------------

### Pipeline parameters

-   **`--input`** (mandatory) The user can provide input in either of the following ways:

    -   A single SRA run accession. eg: **`perl ama.pl --input SRR12548227 --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa`**

    -   A list of run accessions in a text file (1 run accession per line). eg: **`perl ama.pl --input example/list.txt --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa`**

    -   The SRA runInfo exported directly from the NCBI-SRA web portal. Goto the [SRA homepage](https://www.ncbi.nlm.nih.gov/sra "Home - NCBI - SRA") and search for the desired keyword. Export the `SraRunInfo.csv` by clicking 'Send to' =\> File =\> RunInfo). eg: **`perl ama.pl --input example/SraRunInfo.csv --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa`**

-   **`--sequences`** (mandatory) The user should provide a fasta file containing the query sequences.

-   **`--output`** (optional) The output directory to store the results. By default, the output will be stored into the **`results/`** directory of the package. eg: **`perl ama.pl --input example/SraRunInfo.csv --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa --output /src/main/test/`**

-   **`--mode`** (optional) Choose one of the three modes to run the pipeline.

    -   The **`screen`** is the default mode which will only download a fraction of the data-set per SRA-run accession and analyse the file as per the given configuration.

    -   The **`full`** mode will execute the pipeline by downloading the complete fastq file per SRA-run accession.

    -   The **`both`** option searches for samples using a fraction of the data that meet the minimum alignment cutoff from either 'bowtie2' or 'blastn', and then automatically performs alignment by downloading the entire fastq file. eg: **`perl ama.pl --input example/SraRunInfo.csv --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa --output /src/main/test/ --mode screen`**

        > *Note:* There is a supporting **`summary`** mode, that will generate a unified alignment summary by examining the output files created by either screen-mode or full-mode. The summary mode should only be used when the user needs to recreate the summary stats from the pre-existing results. The user must enter **`–mode summary`** along with the previously used command parameters to re-generate the summary.

    -   **`--config`** (optional) Pipeline configuration. By default it will use the **`conf.txt`** generated by the setup script. eg: **`perl ama.pl --input example/SraRunInfo.csv --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa --output /src/main/test/ --mode screen --config conf.txt`**

------------------------------------------------------------------------

### Troubleshooting {#troubleshooting}

-   Errors related to mamba/conda environment:

    Since `mamba` is a drop-in replacement and uses the same commands and configuration options as **conda**, it's possible to swap almost all commands between **conda** & **mamba**.

    Use **`conda list`** command to verify whether the packages mentioned in the `requirements.yaml` are successfully installed into your environment.

    > *Note:* The `requirements.yaml` provided in this package was exported from `mamba 0.25.0` installation running on `Ubuntu 20.04.4 LTS`.

    In case of any missing tool/ conflicting dependencies in the environment, the user can try using **`conda search <tool name>`** or `mamba repoquery search` command to find the supported version of the tool and then manually install it by typing **`conda install <tool name>`** or `mamba install <tool name>` inside the environment. Please refer the official [troubleshooting guide](https://conda.io/projects/conda/en/latest/user-guide/troubleshooting.html "User guide » Troubleshooting") for further help.

    > *Note:* On macOS and Linux, the supported tools and their dependencies aren't always the same. Even when all of the requirements are completely aligned, the set of available versions isn't necessarily the same. User may try setting up the environment using any of the supplementary `requirements-*.txt` provided in the `src/main/resources/` directory.

-   Error installing Perl modules:

    Users must ensure that they have write permission to the `/Users/\*/.cpan/` or similar directory, and the CPAN is properly configured.

    > *Note about MAKE*: 'make' is an essential tool for building Perl modules. Please make sure that you have 'make' installed in your system. The setup script provided in this package utilizes 'cpan' to build the required Perl modules automatically.

------------------------------------------------------------------------

### List of Perl modules and tools incorporated in the pipeline {#list-of-perl-modules-and-tools-incorporated-in-the-pipeline}

-   Perl modules:

    -   Config::Simple
    -   Parallel::ForkManager
    -   Log::Log4perl
    -   Getopt::Long
    -   Text::CSV
    -   Text::Unidecode

-   Tools:

    -   [NCBI EDirect utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
    -   [NCBI SRA Toolkit](https://www.ncbi.nlm.nih.gov/home/tools/)
    -   [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
    -   [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
    -   [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
    -   [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    -   [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    -   [Samtools](http://www.htslib.org/download/)
    -   [Kraken2](https://ccb.jhu.edu/software/kraken2/)

------------------------------------------------------------------------
