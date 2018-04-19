FROM matsengrp/cpp

# ----------------------------------------------------------------------------------------
RUN pip install numpy
RUN pip install pandas
RUN pip install biopython

COPY . /PyMotifFinder
WORKDIR /PyMotifFinder
RUN pip install .
CMD python test/test_motif_finder.py && python test/test_process_partis.py