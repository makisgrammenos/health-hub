# FROM node:22.11-alpine

# # Optional: Fix for watching issues (Alpine is minimal)
# RUN apk add --no-cache bash git

# WORKDIR /app

# COPY package.json ./
# RUN npm install --legacy-peer-deps

# COPY . .

# ENV CHOKIDAR_USEPOLLING=true
# EXPOSE 3000

# CMD ["npm", "run", "dev"]
# healthhub-next/Dockerfile
FROM node:22-bullseye-slim

# Optional but handy
ENV NODE_ENV=development \
    NEXT_TELEMETRY_DISABLED=1 \
    NPM_CONFIG_UPDATE_NOTIFIER=false

# OS deps you might need for git deps / fonts
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Leverage Docker cache properly
COPY package.json  ./
# RUN npm ci --no-audit --no-fund
RUN npm i --legacy-peer-deps
# Now copy the rest
COPY . .

EXPOSE 3000
CMD ["npm", "run", "dev"]
