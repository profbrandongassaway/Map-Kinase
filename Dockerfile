# syntax=docker/dockerfile:1

############################
# Runtime stage
############################
FROM python:3.12-slim AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    M5_HOST=0.0.0.0 \
    M5_PORT=8080

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends netcat-openbsd \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --shell /bin/bash appuser

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=appuser:appuser . /app

RUN mkdir -p /app/cache \
    && chgrp -R 0 /app \
    && chmod -R g=u /app

USER appuser

EXPOSE 8080


CMD ["python", "-m", "uvicorn", "MapKinase_WebApp.m5_main_ui:app", \
     "--host", "0.0.0.0", \
     "--port", "8080", \
     "--log-level", "info"]