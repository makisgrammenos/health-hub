FROM node:22.11-alpine

WORKDIR /app

COPY package.json package-lock.json ./
RUN npm install --legacy-peer-deps

COPY . .

EXPOSE 3000

# Enable file change polling for hot reload
ENV CHOKIDAR_USEPOLLING=true

CMD ["npm", "run", "dev"]
