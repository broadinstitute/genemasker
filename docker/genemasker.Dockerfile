FROM python:3.9.21-slim-bookworm

# all main work in work
WORKDIR /work

# install modules
#pip3 install numpy pandas scipy scikit-learn psutil "dask[dataframe]"
RUN pip3 install --upgrade pip && \
    pip3 install numpy==2.0.2 pandas==2.2.3 scikit-learn==1.6.1 scipy==1.13.1 psutil==6.1.1 dask[dataframe]==2024.8.0

# open /work directory permissions
RUN chmod -R 777 /work
