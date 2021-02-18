ARG BASEIMAGE=ontresearch/base-workflow-image:v0.1.0
FROM $BASEIMAGE

# Minimal install for example purposes
COPY environment.yaml $HOME/environment.yaml 
RUN \
    awk '{if(index($0, "- pip")){exit 0}else{print}}' $HOME/environment.yaml > $HOME/environment.fixed.yaml \
    && . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install --file $HOME/environment.fixed.yaml \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && pip install "aplanat>=0.2.8"

USER $WF_UID
WORKDIR $HOME
