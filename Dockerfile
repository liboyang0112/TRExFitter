FROM atlasamglab/stats-base:root6.23.01

SHELL [ "/bin/bash", "-c" ]

COPY . /code/src/TRExFitter

# we need to remove the build directory as it is leftover from the compile step of the CI
RUN cd /code/src/TRExFitter && \
    rm -rf build && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j4 && \
    cd ..

# add to path to allow running via "trex-fitter"
RUN echo "export PATH=/code/src/TRExFitter/build/bin:${PATH}" >> /home/.bashrc

WORKDIR /home/data
ENV HOME /home

ENTRYPOINT ["/bin/bash", "-l", "-c"]
CMD ["/bin/bash"]
