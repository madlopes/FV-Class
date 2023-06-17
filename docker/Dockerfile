FROM python:3.9.9

# Installing python 2.7
WORKDIR /usr/src
RUN apt-get update \
    && apt-get install build-essential checkinstall \
    && wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz \
    && tar xzf Python-2.7.18.tgz \
    && cd Python-2.7.18 \
    && ./configure --enable-optimizations \
    && make altinstall \
    && wget https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2.7 get-pip.py
    
    
WORKDIR /usr/bin/
COPY ./docker/artifacts  artifacts/
RUN chmod +x artifacts/probe/probe.2.16r72.130520 artifacts/reduce/reduce.3.23.130521 \
    && ln -s artifacts/probe/probe.2.16r72.130520 probe \
    && ln -s artifacts/reduce/reduce.3.23.130521 reduce
    
WORKDIR /main
COPY ./docker/requirements  requirements/
RUN pip2 install -r requirements/python2_requirements.txt \
    && pip3 install -r requirements/python3_requirements.txt

RUN wget -c https://rinalyzer.de/rinerator/RINerator_V0.5.1.tar.gz -O - | tar -xz
    # && tar xzf RINerator_V0.5.1.tar.gz --strip-components=1 \
    # && rm RINerator_V0.5.1.tar.gz

COPY ./data data/
COPY ./src src/
