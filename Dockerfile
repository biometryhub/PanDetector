# Dockerfile to automate setup of the PanDetector multiple genome
# alignment tool into a working out-of-box Docker container.
#
# Assuming that you are running this Dockerfile from the GitHub
# working directory and that your genome data files are stored
# under the example/ subdirectory, build and run this container
# using Docker/Podman with:
#   docker build -t pandetector .
#   docker run -it -v $(pwd)/example:/example pandetector
#
# The genome aligner can then be run in the expected way in the 
# container shell:
#   ./pipeline_Minimap.sh -g example/quick_genomes.txt -o example/output -d 20
#
# Code author: Russell A. Edson, AAGI-AU/Biometry Hub
# Date last modified: 23/09/2026
FROM docker.io/ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && apt-get -y install \
  build-essential \
  git \
  wget \
  openjdk-17-jdk \
  zlib1g-dev

# Download + install Minimap2-2.31 release
RUN wget "https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31.tar.bz2" \
  && tar -xjf minimap2-2.31.tar.bz2 \
  && cd minimap2-2.31 \
  && make \
  && cp minimap2 misc/paftools.js /usr/bin/ \
  && cd ..

# Download + install k8-1.2 release
RUN wget "https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2" \
  && tar -xjf k8-1.2.tar.bz2 \
  && cp k8-1.2/k8-x86_64-Linux /usr/bin/k8

# Download + install GSAlign (latest version from git repository)
RUN git clone "https://github.com/hsinnan75/GSAlign" \
  && cd GSAlign \
  && make \
  && cp bin/bwt_index bin/GSAlign /usr/bin/ \
  && cd ..

# Pull PanDetector and install
RUN mkdir -p /PanDetector
COPY . /PanDetector/
RUN cp PanDetector/MFbio.jar PanDetector/pipeline_Minimap.sh \
  PanDetector/pipeline_GSAlign.sh / \
  && chmod u+x pipeline_Minimap.sh pipeline_GSAlign.sh

# Helpful message for container startup
RUN echo 'echo -e "Run the PanDetector Multiple Genome Alignment tool using:\n"\
  "  ./pipeline_Minimap.sh -g example/quick_genomes.txt -o example/output"\
  "-d 20 -n 4\nfor example. Refer to the usage guide available at"\
  "the GitHub page (in README.md and Manual.md) for the options and input formats."\
  "\nNote: make sure to write the output files to a directory under the"\
  "example/ bind-mount, so that you can access them from outside of the"\
  "Docker container afterward!\n"' > /etc/profile.d/welcome.sh

CMD [ "/bin/bash", "-l" ]
