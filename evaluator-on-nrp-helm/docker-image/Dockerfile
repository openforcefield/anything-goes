FROM condaforge/mambaforge:23.3.1-1


SHELL ["/bin/bash", "-c"]

ENV PATH /opt/conda/bin:$PATH
ENV PYTHON_VERSION=${python}
ENV DASK_VERSION=${release}

RUN mamba install -y -c conda-forge \
    "mamba>=0.27.0" \
    python=3.10 \
    nomkl \
    cmake \
    dask=2024.2.1 \
    cachey \
    streamz \
    openff-evaluator==0.4.9 \
    rdkit       \
    ambertools   \
    && mamba clean -tipy \
    && find /opt/conda/ -type f,l -name '*.a' -delete \
    && find /opt/conda/ -type f,l -name '*.pyc' -delete \
    && find /opt/conda/ -type f,l -name '*.js.map' -delete \
    && find /opt/conda/lib/python*/site-packages/bokeh/server/static -type f,l -name '*.js' -not -name '*.min.js' -delete \
    && rm -rf /opt/conda/pkgs

RUN mkdir -p ~/.config/dask
COPY distributed.yaml ~/.config/dask/distributed.yaml
COPY server.py /usr/bin/server.py
COPY run-server.sh /usr/bin/run-server.sh

RUN mkdir /opt/app

ENTRYPOINT ["tini", "-g", "--", "/usr/bin/run-server.sh"]
