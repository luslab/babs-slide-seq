#FROM nfcore/base:2.1
#FROM centos:centos7
#FROM nvidia/cuda:11.0.3-devel-centos7
FROM nvidia/cuda:11.0.3-cudnn8-devel-centos7

#RUN yum -y update

RUN yum -y install gcc gcc-c++ make kernel-devel lshw
RUN yum -y install cuda-nsight-systems-10-1-10.1.168-1.x86_64

ENV LD_LIBRARY_PATH /usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64
ENV LD_LIBRARY_PATH /home/hamming/opencl/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH /usr/local/cuda-10.1/NsightSystems-2019.3/Host-x86_64:$LD_LIBRARY_PATH

# Copy files over
ADD hamming /home/hamming

WORKDIR /home/hamming
RUN make hamming

# RUN apt-get update \
#  && apt-get install -y --no-install-recommends \
#  procps=2:3.3.15-2 \
#  libz-dev \
#  build-essential \
#  && apt-get clean \
#  && rm -rf /var/lib/apt/lists/*

# RUN mkdir -p /.config/matplotlib && chmod 777 /.config/matplotlib
# RUN mkdir -p /.cache/matplotlib && chmod 777 /.cache/matplotlib

# # Update conda
# RUN conda update -n base -c defaults conda

# # Install conda packages
# COPY environment.yml /
# RUN conda env create -f /environment.yml && conda clean -a
# ENV PATH /opt/conda/envs/nfcore-module-tobias/bin:$PATH

# Install dev fork
# WORKDIR /home
# RUN git clone https://github.com/chris-cheshire/TOBIAS.git
# WORKDIR /home/TOBIAS
# RUN /opt/conda/envs/nfcore-module-tobias/bin/pip install .
# WORKDIR /home
# RUN conda install -n nfcore-module-tobias -c conda-forge filprofiler