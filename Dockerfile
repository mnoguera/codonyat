FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

COPY pyproject.toml poetry.lock* /app/

RUN python -m pip install --upgrade pip setuptools wheel
RUN pip install pyproject.toml poetry.lock* --no-cache-dir

CMD ["codonyat"]
