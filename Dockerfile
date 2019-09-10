FROM python:3.7

# ----------------------------------------------------------------------------------------
RUN pip install biopython
RUN pip install pandas

COPY . /PyMotifFinder
WORKDIR /PyMotifFinder
RUN pip install .
CMD python test/test_motif_finder.py && python test/test_process_partis.py
