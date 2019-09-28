FROM placeholder-will-be-replaced-with-CI-variable

## ensure locale is set during build
ENV LANG C.UTF-8

COPY . /TRexFitter/source/TRexFitter/
WORKDIR /TRexFitter

USER root

RUN cd /TRexFitter && \
    cp source/TRexFitter/util/asetup-CMakeLists.txt source/CMakeLists.txt && \
    mkdir /TRexFitter/build && \
    cd build && \
    source /release_setup.sh && \
    cmake ../source/. && \
    make

CMD source /release_setup.sh && source  /TRexFitter/build/x86_64-centos7-gcc8-opt/setup.sh && /bin/bash
