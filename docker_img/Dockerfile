FROM ubuntu:16.04
LABEL mantainer="Rodrigo Theodoro Rocha"
LABEL git_mantainer="github.com/rodtheo"

RUN apt-get update && apt-get install -y r-base libpath-tiny-perl git wget g++ python2.7 make libz-dev

RUN cpan -i File::Basename
RUN cpan -i File::Copy
RUN cpan -i File::Spec
RUN cpan -i File::Spec::Link
RUN cpan -i Getopt::Long
RUN cpan -i List::Util

WORKDIR home

RUN wget ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.18.tar.gz
RUN tar xvzf Reapr_1.0.18.tar.gz
WORKDIR Reapr_1.0.18
WORKDIR /home/Reapr_1.0.18/src

RUN cp task_preprocess.pl task_preprocess2.pl
RUN sed '343s|^.*$|system_call("Rscript $r_script >> " . $r_script . "out");|' task_preprocess2.pl > task_preprocess.pl

RUN cp task_fcdrate.cpp task_fcdrate2.cpp
RUN sed '208s|^.*$|//systemCall("Rscript --no-save " + outfile + " > " + outfile + "out");|'  task_fcdrate2.cpp > task_fcdrate.cpp

RUN cp histogram.cpp histogram2.cpp
RUN sed '381s|^.*$|//systemCall("Rscript " + outfile_R + " > " + outfile_R + "out");|' histogram2.cpp > histogram.cpp

WORKDIR /home/Reapr_1.0.18/
RUN ./install.sh
ENV PATH="/home/Reapr_1.0.18/:${PATH}"

WORKDIR /home
RUN git clone https://github.com/sc932/ALE.git
WORKDIR ALE/src
RUN make
ENV PATH="/home/ALE/src:${PATH}"

RUN apt-get -qy autoremove
WORKDIR /home

#RUN export PATH="/home/Reapr_1.0.18:$PATH
