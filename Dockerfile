# Dockerfile
FROM ghcr.io/prefix-dev/pixi:0.23.0 AS build

COPY . /app
WORKDIR /app
RUN pixi install
RUN pixi run install
RUN pixi shell-hook > /shell-hook
RUN chmod +x /shell-hook

WORKDIR /app
ENTRYPOINT [ "pixi", "run" ]

