# Base image is pinned by tag and digest; the image is cached locally and the
# build succeeds with --network=none.
FROM python:3.11.17-slim@sha256:0dd364ba7e10242f07755449e3a3d0e35f9efd987952737b90def6709ab0c5ce

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Agent identity is set by [agent].user in task.toml. Do not switch USER here:
# Harbor's keep-alive process must remain owned by the image default user.
RUN useradd --create-home --uid 1000 runner

WORKDIR /app
COPY app/ /app/
RUN chown -R runner:runner /app
