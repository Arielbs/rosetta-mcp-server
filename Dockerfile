FROM node:22-slim

# Install Python
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv && rm -rf /var/lib/apt/lists/*

# Create venv and install biotite
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir biotite

# Set Python path for the MCP server
ENV PYTHON_BIN=/opt/venv/bin/python

WORKDIR /app
COPY package.json rosetta_mcp_wrapper.js rosetta_mcp_http.js rosetta_mcp_server.py ./
COPY docs_cache/ ./docs_cache/

EXPOSE 8080
CMD ["node", "rosetta_mcp_http.js"]
