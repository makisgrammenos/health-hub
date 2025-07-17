FROM node:22.11-alpine

# Optional: Fix for watching issues (Alpine is minimal)
RUN apk add --no-cache bash git

WORKDIR /app

COPY package.json ./
RUN npm install --legacy-peer-deps

COPY . .

ENV CHOKIDAR_USEPOLLING=true
EXPOSE 3000

CMD ["npm", "run", "dev"]
