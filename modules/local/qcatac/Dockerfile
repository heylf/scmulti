# Use an official Python runtime as a parent image
FROM python:3.10

# Install packages
RUN pip install snapatac2
RUN pip install nbformat
RUN pip install plotly
RUN pip install muon

# Set the environment variable NUMBA_CACHE_DIR
ENV NUMBA_CACHE_DIR '/tmp'

# Set the working directory to /app
WORKDIR /home/snapatac2

# Copy the current directory contents into the container at /app
COPY . /home/snapatac2

# Run app.py when the container launches
CMD ["python3"]
