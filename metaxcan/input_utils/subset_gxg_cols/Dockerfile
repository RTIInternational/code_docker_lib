FROM ubuntu:18.04

# Add scripts to make it run
ADD subset_gxg_cols.sh /
ADD run.sh /

# Add to /opt folder
RUN mkdir /opt/code_docker_lib

# Move scripts to
RUN mv subset_gxg_cols.sh /opt/code_docker_lib && \
    mv run.sh /opt/code_docker_lib

# Change permissions to make things exectuable
RUN chmod -R u+x /opt/code_docker_lib/

# Set wrkdir
RUN mkdir /data
WORKDIR /data

# Create entrypoint that runs the run.sh wrapper script
ENTRYPOINT ["/bin/bash", "/opt/code_docker_lib/run.sh"]