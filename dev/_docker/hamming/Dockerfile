FROM nvidia/cuda:11.0.3-runtime-centos7

# Install packages
RUN yum -y install gcc gcc-c++ make kernel-devel lshw
RUN yum -y install cuda-nsight-systems-10-1-10.1.168-1.x86_64

# Set lib paths
ENV LD_LIBRARY_PATH /usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64
ENV LD_LIBRARY_PATH /home/hamming/opencl/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH /usr/local/cuda-10.1/NsightSystems-2019.3/Host-x86_64:$LD_LIBRARY_PATH

# Copy files over
ADD hamming /home/hamming

# Build files
WORKDIR /home/hamming
RUN make hamming