# PanDetector Multiple Genome Aligner
<!-- badges -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub License](https://img.shields.io/github/license/biometryhub/PanDetector)
![Static Badge](https://img.shields.io/badge/version-1.0.0-80b6ff)

PanDetector is a new fast and flexible program that is able to identify the pan-genome sequence of larger and more evolutionary diverse genomes. 

- [Quick start](#qstart)
- [Quick start (using Docker)](#dockerqstart)
- [Usage](#usage)

## <a name="qstart"></a>Quick start
Installation and configuration of PanDetector on Linux-based operating systems proceeds as follows.

#### Step 1. Configure your `$PATH` for PanDetector binary dependencies
PanDetector depends on the [Minimap2](https://github.com/lh3/minimap2) versatile pairwise aligner (and its related `paftools.js` utility), as well as the [K8 Javascript shell](https://github.com/attractivechaos/k8). The easiest way is to install these to a prepared folder on the system `$PATH` for them, so that they are always available when PanDetector runs:
```bash
mkdir -p $HOME/bin
echo "export PATH=$HOME/bin:${PATH}" >> $HOME/.bashrc && source $HOME/.bashrc
```

#### Step 2. Download and install Minimap2 (v2.26)
Grab the v2.26 release of Minimap2 from its GitHub repository [here](https://github.com/lh3/minimap2/releases/tag/v2.26). Alternatively, copy-paste the below commands to automatically download, compile and configure Minimap2. (Note: this compilation requires compiler tools and the zlib development headers to be installed: on Ubuntu 22.04, you can easily install these compilation dependencies with `sudo apt-get -y install build-essential zlib1g-dev`. You might need to run `sudo apt-get update` before installing build-essential and zlib1g-dev)
```bash
wget "https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26.tar.bz2"
tar -xjf minimap2-2.26.tar.bz2
cd minimap2-2.26 && make
cp minimap2/misc/paftools.js $HOME/bin/
cd ..
```

#### Step 3. Download and install K8 (v1.0)
Get the v1.0 release of the K8 Javascript shell from its GitHub repository [here](https://github.com/attractivechaos/k8/releases/tag/v1.0). Alternatively, execute the following commands to automatically download and configure the precompiled K8 binary:
```bash
wget "https://github.com/attractivechaos/k8/releases/download/v1.0/k8-1.0.tar.bz2"
tar -xjf k8-1.0.tar.bz2
cp k8-1.0/k8-x86_64-Linux $HOME/bin/k8
```

#### Step 4. Install a Java runtime/development kit
[OpenJDK-11](https://openjdk.org/projects/jdk/11/) (or later versions) have been confirmed to work well with PanDetector. For most Linux systems, these are easily installed via the package manager. E.g., to install OpenJDK-11 (the default JDK) on Ubuntu 22.04:
```bash
sudo apt-get -y install openjdk-11-jdk  # or default-jdk
```

#### Step 5. Download PanDetector and run an example pipeline
Finally, pull this GitHub repository to download the PanDetector tool, and run a test case on the provided example set of genomes to confirm that the tool is working correctly.
```bash
git clone https://github.com/mfruzan/PanDetector.git
cd PanDetector
chmod +x pipeline_Minimap.sh

./pipeline_Minimap.sh -g example/quick_genomes.txt -o example_out -d 20 -n 16
```

## <a name="dockerqstart"></a>Quick start (using Docker)
Alternatively, easily set up PanDetector in a Docker container using the provided Dockerfile, which completely automates the installation. For information about setting up Docker on Windows/Mac/Linux and using containers, see [docs.docker.com](https://docs.docker.com/).
```bash
git clone https://github.com/mfruzan/PanDetector.git
cd PanDetector
sudo docker build -t Pandetector .
sudo docker run -it -v $(pwd)/example:/example pandetector
```
In the interactive shell for the container, you can immediately run the Multiple Genome Aligner tool:
```bash
./pipeline_Minimap.sh -g example/quick_genomes.txt -o example/output -d 20 -n 16
```

## <a name="usage"></a>Usage
Use the PanDetector multiple alignment tool (with the Minimap2 pipeline) as follows:
```bash
./pipeline_Minimap.sh -g <genome_list> -o <out_dir> -d <divergence> -n <ncpus> -m <minlength> -c <chromosome>
```

The main input file for PanDetector is the `<genome_list>` text file, consisting of lines of genomes:
```bash
Alg130	example/Alg130.fna
DW5	example/DW5.fna
M4	example/M4.fna
```
Each line contains an alias name (e.g., Alg130, DW5), followed by a space/Tab, then followed by the filepath to the FASTA file for that genome. In this example, Alg130 is the query genome, and the rest of the genomes become the subjects. This text file is passed to `./pipeline_Minimap.sh` using the `-g` flag.

The `-o` argument specifies the output directory. Note that this
directory will be created if it does not already exist. Inside the
output directory, PanDetector generates a series of intermediate output files during the alignment process, but the main program outputs are:

- `core_msa.maf.gz`: This is a gzipped version of a standard MAF file, each entry of core contains subject lines related to all the genomes. Coordinates and strandness of entries are in respect to the original genome FASTA file. 
- `accessory_msa.maf.gz`: This is a gzipped version of a standard MAF file, each entry of accessory contains subject lines related to some (but not all) of the genomes. Coordinates and strandness of entries are in respect to the original genome FASTA file. 

The `-d` argument is the expected divergence level, and can be any integer between 1 and 40.

Other arguments to PanDetector are optional, and allow fine-tuning of the program configuration:

- `-n` is the number of cores/CPUs to use for the program execution (default is 4 cores).
- `-m` is the minimum alignment length, in bp (the default is 200bp).
- `-c` toggles chromosome number matching (1: enabled, 0:disabled, default is 0). Note that if chromosome number matching is enabled, PanDetector considers a contig name to start with a chromosome number, such as '2B' or 'chr14' or simply '14', followed by a space (or the characters '_','-','!'). If this pattern does not exist in contig names, no error will be raised, but chromosome checking will be skipped for current genome.

The [PanDetector Manual](https://github.com/biometryhub/PanDetector/blob/master/Manual.md) explains program usage in detail, and lists further analysis examples.
