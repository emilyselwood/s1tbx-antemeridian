FROM continuumio/miniconda3:4.7.10

LABEL maintainer="Emily Selwood"

RUN apt-get update \
    && apt-get install -y build-essential libgdal-dev \
    && wget --quiet http://step.esa.int/downloads/7.0/installers/esa-snap_sentinel_unix_7_0.sh \
    && /bin/sh ./esa-snap_sentinel_unix_7_0.sh -q \
    && rm ./esa-snap_sentinel_unix_7_0.sh \
    && /opt/snap/bin/snap --nosplash --nogui --modules --update-all \
    && pip install --no-cache-dir \
    asynchronousfilereader \
    xmltodict \
    scipy \
    numpy \
    pyproj==1.9.6 \
    GDAL==2.4.0 \
    rasterio \
    matplotlib \
    redis

ADD . /app/
WORKDIR /app/

COPY snap/bin/gpt.vmoptions /opt/snap/bin/
COPY snap/etc/snap.auxdata.properties /opt/snap/etc/

CMD python run.py S1A_IW_GRDH_1SDV_20191001T173212_20191001T173241_029268_03535A_5270.zip --out_path /data/ && tail -f /var/log/*.log
