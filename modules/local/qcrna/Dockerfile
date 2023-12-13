# Use an official Python runtime as a parent image
FROM python:3.10

# Install packages
RUN pip install muon
RUN pip install nbformat
RUN pip install pysam
RUN pip install plotly
RUN pip install -U numba
RUN pip install -U scanpy

# Set the environment variable NUMBA_CACHE_DIR
ENV NUMBA_CACHE_DIR '/tmp'

# Set the working directory to /app
WORKDIR /home/muon

# Copy the current directory contents into the container at /app
COPY . /home/muon

# Run app.py when the container launches
CMD ["python3"]
