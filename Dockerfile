# syntax=docker/dockerfile:1

############################
# Builder stage
############################
FROM registry.access.redhat.com/hi/python:latest-builder AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app
USER root

COPY requirements.txt /app/

RUN python -m pip install --prefix=/install -r requirements.txt

COPY . /app

# Directories that OpenShift can write to
RUN mkdir -p /app/cache \
    && chgrp -R 0 /app \
    && chmod -R g=u /app

############################
# Runtime stage
############################
FROM registry.access.redhat.com/hi/python:latest AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    M5_HOST=0.0.0.0 \
    M5_PORT=8080

WORKDIR /app
USER root

# Install netcat, matching the dependency from the original image.
RUN microdnf install -y nmap-ncat \
    && microdnf clean all

# Copy installed Python packages and application source from builder.
COPY --from=builder /install /usr/local
COPY --from=builder /app /app

# Support OpenShift's arbitrary, non-root runtime UID.
RUN chgrp -R 0 /app \
    && chmod -R g=u /app

USER 1001

EXPOSE 8080

CMD ["python", "-m", "uvicorn", \
    "MapKinase_WebApp.m5_main_ui:app", \
    "--host", "0.0.0.0", \
    "--port", "8080", \
    "--log-level", "info", \
    "--timeout-graceful-shutdown", "30", \
    "--no-access-log"]