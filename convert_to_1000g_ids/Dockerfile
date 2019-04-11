FROM ubuntu:18.04

# Add scripts to make it run
ADD convert_to_1000g_ids.pl /
ADD run.sh /

# Add to /opt folder
RUN mkdir /opt/code_docker_lib
RUN mv convert_to_1000g_ids.pl /opt/code_docker_lib

# Change permissions to make things exectuable
RUN chmod u+x /opt/code_docker_lib/convert_to_1000g_ids.pl

# Set wrkdir
RUN mkdir /data
WORKDIR /data

# Create entrypoint that runs the run.sh wrapper script
# ENTRYPOINT ["/bin/bash", "/opt/code_docker_lib/run.sh"]