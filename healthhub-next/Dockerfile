# Use official Node.js image
FROM node:18

# Set the working directory
WORKDIR /app

# Copy package.json and install dependencies
COPY package.json ./
RUN npm install

# Copy the source code
COPY ./ ./

# Expose the Next.js port
EXPOSE 3000

# Start the Next.js server in development mode
CMD ["npm", "run", "dev"]
