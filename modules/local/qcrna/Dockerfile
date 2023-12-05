# Use an official Python runtime as a parent image
FROM python:3.10

# Install Muon
RUN pip install muon

# Set the working directory to /app
WORKDIR /home/muon

# Copy the current directory contents into the container at /app
COPY . /home/muon

# Run app.py when the container launches
CMD ["python3"]
