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
# FROM python:3.11-slim

# # system deps for opencv
# RUN apt-get update && apt-get install -y --no-install-recommends libgl1 && rm -rf /var/lib/apt/lists/*

# WORKDIR /app
# COPY req.txt .
# RUN pip install  --upgrade pip && pip install -r req.txt

# COPY . .
# EXPOSE 8000
# CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000","--reload"]


# syntax=docker/dockerfile:1.7
FROM python:3.11-slim AS base
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Minimal system deps (opencv needs libgl1)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install deps first (cacheable layer)
# Prefer a single pinned file: requirements.txt (rename req.txt -> requirements.txt)
COPY requirements.txt .
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip && \
    pip install -r requirements.txt

# Copy only code, NOT models/data/uploads
COPY app.py ./app.py
COPY routers/ ./routers/
COPY utils/ ./utils/

EXPOSE 8000
# NOTE: use --reload only in dev (compose). Keep prod command clean.
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]

