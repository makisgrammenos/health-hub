# FROM python:3.11-slim

# WORKDIR /app
# # Build deps (for annoy and other C-extensions) + a couple of common runtime libs
# RUN apt-get update && apt-get install -y --no-install-recommends \
#     build-essential gcc g++ python3-dev \
#     libglib2.0-0 libgl1 git \
#  && rm -rf /var/lib/apt/lists/*
# COPY req.txt .
# RUN pip install  -r req.txt

# COPY . .

# CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7001","--reload"]

# syntax=docker/dockerfile:1.7
# syntax=docker/dockerfile:1.7
# syntax=docker/dockerfile:1.7
FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ python3-dev \
    libglib2.0-0 libgl1 git \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip && \
    pip install -r requirements.txt

# Copy only code; datasets will be mounted
COPY *.py ./

EXPOSE 7001
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7001"]
