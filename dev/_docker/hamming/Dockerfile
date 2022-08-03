# Taken from https://gitlab.com/nvidia/container-images/cuda/-/tree/master/dist/11.4.0/centos7
################## BASE ##################
FROM centos:7 as base

FROM base as base-amd64

ENV NVARCH x86_64
ENV NVIDIA_REQUIRE_CUDA "cuda>=11.4 brand=tesla,driver>=418,driver<419 brand=tesla,driver>=450,driver<451"
ENV NV_CUDA_CUDART_VERSION 11.4.43-1

COPY cuda.repo-x86_64 /etc/yum.repos.d/cuda.repo

FROM base-${TARGETARCH}

ARG TARGETARCH

LABEL maintainer "NVIDIA CORPORATION <sw-cuda-installer@nvidia.com>"

RUN NVIDIA_GPGKEY_SUM=d0664fbbdb8c32356d45de36c5984617217b2d0bef41b93ccecd326ba3b80c87 && \
    curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/rhel7/${NVARCH}/D42D0685.pub | sed '/^Version/d' > /etc/pki/rpm-gpg/RPM-GPG-KEY-NVIDIA && \
    echo "$NVIDIA_GPGKEY_SUM  /etc/pki/rpm-gpg/RPM-GPG-KEY-NVIDIA" | sha256sum -c --strict -

ENV CUDA_VERSION 11.4.0

# For libraries in the cuda-compat-* package: https://docs.nvidia.com/cuda/eula/index.html#attachment-a
RUN yum upgrade -y && yum install -y \
    cuda-cudart-11-4-${NV_CUDA_CUDART_VERSION} \
    cuda-compat-11-4 \
    && ln -s cuda-11.4 /usr/local/cuda \
    && yum clean all \
    && rm -rf /var/cache/yum/*

# nvidia-docker 1.0
RUN echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf

ENV PATH /usr/local/nvidia/bin:/usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/local/nvidia/lib:/usr/local/nvidia/lib64

#COPY NGC-DL-CONTAINER-LICENSE /

# nvidia-container-runtime
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility

################## RUN TIME ##################

ENV NV_NVTX_VERSION 11.4.43-1
ENV NV_LIBNPP_VERSION 11.4.0.33-1
ENV NV_LIBNPP_PACKAGE libnpp-11-4-${NV_LIBNPP_VERSION}
ENV NV_LIBCUBLAS_VERSION 11.5.2.43-1
ENV NV_LIBNCCL_PACKAGE_NAME libnccl
ENV NV_LIBNCCL_PACKAGE_VERSION 2.11.4-1
ENV NV_LIBNCCL_VERSION 2.11.4
ENV NCCL_VERSION 2.11.4
ENV NV_LIBNCCL_PACKAGE ${NV_LIBNCCL_PACKAGE_NAME}-${NV_LIBNCCL_PACKAGE_VERSION}+cuda11.4

RUN yum install -y \
    cuda-libraries-11-4-${NV_CUDA_LIB_VERSION} \
    cuda-nvtx-11-4-${NV_NVTX_VERSION} \
    ${NV_LIBNPP_PACKAGE} \
    libcublas-11-4-${NV_LIBCUBLAS_VERSION} \
    ${NV_LIBNCCL_PACKAGE} \
    && yum clean all \
    && rm -rf /var/cache/yum/*

################## OPEN CL ##################

# Configure NVIDIA/CUDA OpenCL settings
RUN mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

RUN \
    ln -s /usr/local/cuda/lib64/libOpenCL.so /usr/lib64/libOpenCL.so 

# Library Path at BUILD time
ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}
ENV LD_LIBRARY_PATH /usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64:$LD_LIBRARY_PATH

#RUN \
#  ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/local/cuda/lib64/libcuda.so.1 && \
#  ldconfig

ENV CUDA_HOME /usr/local/cuda
ENV OPENCL_LIB /usr/local/cuda/lib64/

#install additional tools and library prerequisites for additional packages
RUN \
  yum install -y epel-release \
  yum install -y opencl-headers

################## BUILD ##################

# Install build tools
RUN yum -y install gcc gcc-c++ make kernel-devel lshw

# Install fixes to /usr/lib64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by hamming)
RUN yum -y install cuda-nsight-systems-10-1-10.1.168-1.x86_64
ENV LD_LIBRARY_PATH /usr/local/cuda-10.1/NsightSystems-2019.3/Host-x86_64:$LD_LIBRARY_PATH

# Copy files over
ADD hamming /home/hamming

# Build files
WORKDIR /home/hamming
RUN make hamming

ENV PATH /home/hamming:${PATH}










# FROM nvidia/cuda:11.0.3-runtime-centos7

# # Install packages
# RUN yum -y install gcc gcc-c++ make kernel-devel lshw
# RUN yum -y install cuda-nsight-systems-10-1-10.1.168-1.x86_64
# RUN yum -y install epel-release
# RUN yum -y install opencl-headers
# RUN yum -y install ocl-icd-devel-2.2.12-1.el7.x86_64

# RUN ln -s /usr/lib64/libOpenCL.so.1 /usr/lib/libOpenCL.so

# # Set lib paths
# ENV LD_LIBRARY_PATH /usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64
# #ENV LD_LIBRARY_PATH /home/hamming/opencl/lib:$LD_LIBRARY_PATH
# ENV LD_LIBRARY_PATH /usr/local/cuda-10.1/NsightSystems-2019.3/Host-x86_64:$LD_LIBRARY_PATH

# # Copy files over
# ADD hamming /home/hamming

# # Build files
# WORKDIR /home/hamming
# RUN make hamming