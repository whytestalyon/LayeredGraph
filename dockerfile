FROM python:3

COPY *.py /LayeredGraph/
COPY *.sh /LayeredGraph/
COPY templates /LayeredGraph/templates/
COPY HPO_graph_data /LayeredGraph/HPO_graph_data/
RUN ls /LayeredGraph/*
WORKDIR /LayeredGraph
RUN sh install.sh

ENTRYPOINT ["python3", "application.py"]