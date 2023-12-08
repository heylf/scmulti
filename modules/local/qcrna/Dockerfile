# Use an official Python runtime as a parent image
FROM python:3.10

# Install packages
RUN pip install muon
RUN pip install nbformat
RUN pip install pysam
RUN pip install plotly

# Set the working directory to /app
WORKDIR /home/muon

# Copy the current directory contents into the container at /app
COPY . /home/muon

# Run app.py when the container launches
CMD ["python3"]
