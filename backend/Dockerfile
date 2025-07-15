# backend/Dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app files
COPY . .

# Run FastAPI app (adjust `app.py` and app object name if needed)
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
