FROM debian:testing
RUN apt-get -y update && apt-get install -y --no-install-recommends \
	muscle \
	vsearch \
	python-biopython \
	python-numpy \
	python-scipy \
	python-numpy \
	python-pandas \
        curl \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
ENV p_MAJOR=0
ENV p_MINOR=4
ENV p_PATCH=8
ENV p_ARCHIVE=pRESTO_v${p_MAJOR}.${p_MINOR}.${p_PATCH}.tar.gz
RUN curl -SL \
	"http://clip.med.yale.edu/presto/downloadHelper.php?dbTable=presto&downloadFile=${p_ARCHIVE}" \
    | tar xzC /usr/local
