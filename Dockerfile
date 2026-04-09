FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

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
    && pip install --no-cache-dir codonyat==1.0.0 multiqc==1.22.2 \
    && rm -rf /var/lib/apt/lists/*

CMD ["bash"]
