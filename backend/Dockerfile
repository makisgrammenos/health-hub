# FROM python:3.11.0-slim

# # Add system packages required for building wheels
# RUN apt-get update && apt-get install -y --no-install-recommends \
#     build-essential gcc g++ python3-dev \
#     libglib2.0-0 libgl1 git \
#  && rm -rf /var/lib/apt/lists/*
# WORKDIR /app

# COPY requirements.txt .
# RUN pip install --upgrade pip setuptools wheel
# RUN pip install -r requirements.txt

# COPY . .

# CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
FROM python:3.11-slim

# system deps for opencv
RUN apt-get update && apt-get install -y --no-install-recommends libgl1 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY req.txt .
RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir -r req.txt

COPY . .
EXPOSE 8000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
