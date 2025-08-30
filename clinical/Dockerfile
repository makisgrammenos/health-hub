# FROM python:3.11-slim

# WORKDIR /app
# RUN apt-get update && apt-get install -y --no-install-recommends \
#     build-essential gcc g++ python3-dev \
#     libglib2.0-0 libgl1 git \
#  && rm -rf /var/lib/apt/lists/*
# COPY req.txt .
# RUN pip install  -r req.txt

# COPY . .

# CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7000","--reload"]


# syntax=docker/dockerfile:1.7
FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1

# RUN apt-get update && apt-get install -y --no-install-recommends \
#     libgl1 && \
#     libgomp1 && \
#     rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y --no-install-recommends libgl1 && \
    if apt-cache show libgomp1 >/dev/null 2>&1; then \
        apt-get install -y --no-install-recommends libgomp1; \
    else \
        # Fallback: gcc pulls libgomp runtime in Debian/Ubuntu
        apt-get install -y --no-install-recommends gcc; \
    fi && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
COPY uci_hf_df.csv .
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip && \
    pip install -r requirements.txt

# Copy only the Python modules you need to run
COPY *.py ./

EXPOSE 7000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7000"]
