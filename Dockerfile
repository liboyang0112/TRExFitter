FROM placeholder-will-be-replaced-with-CI-variable

## ensure locale is set during build
ENV LANG C.UTF-8

COPY . /TRExFitter/source/TRExFitter/
WORKDIR /TRExFitter

USER root

RUN cd /TRExFitter && \
    cp source/TRExFitter/util/asetup-CMakeLists.txt source/CMakeLists.txt && \
    mkdir /TRExFitter/build && \
    cd build && \
    source /release_setup.sh && \
    cmake ../source/. && \
    make

CMD source /release_setup.sh && source  /TRExFitter/build/x86_64-centos7-gcc8-opt/setup.sh && /bin/bash
