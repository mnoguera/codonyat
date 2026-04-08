FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY aa_caller/ aa_caller/

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bowtie2 \
        default-jre-headless \
        fastp \
        fastqc \
        picard-tools \
        samtools \
    && mkdir -p /usr/picard \
    && ln -sf "$(find /usr/share/java -name 'picard*.jar' | head -n 1)" /usr/picard/picard.jar \
    && python -m pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir . multiqc==1.22.2 \
    && rm -rf /var/lib/apt/lists/*

CMD ["bash"]
