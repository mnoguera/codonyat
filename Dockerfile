FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY aa_caller/ aa_caller/

RUN python -m pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir .

ENTRYPOINT ["codonyat"]
