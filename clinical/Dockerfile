FROM python:3.11-slim

WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ python3-dev \
    libglib2.0-0 libgl1 git \
 && rm -rf /var/lib/apt/lists/*
COPY req.txt .
RUN pip install  -r req.txt

COPY . .

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7000"]
