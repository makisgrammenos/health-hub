# Use official Python image
FROM python:3.10.15-slim

# Install system dependencies for hot reloading
RUN apt-get update && apt-get install -y build-essential libpq-dev

# Set the working directory
WORKDIR /app

# Copy requirements and install dependencies
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir watchdog  # Required for hot reload

# Copy the entire app
COPY . /app

# Expose the FastAPI port
EXPOSE 8000

# Run the FastAPI server with hot reload
CMD ["fastapi", "run", "--reload", "--host", "0.0.0.0", "--port", "8000"]
