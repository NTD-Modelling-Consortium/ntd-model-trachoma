FROM python 
LABEL maintainer="Mike Irvine" 
LABEL version="1.0" 
LABEL description="docker image for Python Trachoma project"  

COPY /.devcontainer/requirements.txt /tmp/
RUN pip install --requirement /tmp/requirements.txt
COPY . /tmp/
RUN pip install -e /tmp

EXPOSE 8888  
CMD ["/bin/bash"]